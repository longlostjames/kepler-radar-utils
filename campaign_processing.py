#!/usr/bin/env python
# coding: utf-8

"""
campaign_processing.py

Campaign-specific processing functions for Kepler MIRA-35 radar data

This module provides processing functions for different campaigns:
- WOEST (WESCON Observing the Evolving Structures of Turbulence)
- CCREST (Characterising CiRrus and icE cloud acrosS the specTrum - Microwave)
- COBALT (Contrail OBservations And Lifecycle Tracking)
- KASBEX (Ka- and S-Band EXperiment)

Author: Chris Walden, UK Research & Innovation and
        National Centre for Atmospheric Science
Last modified: 10-08-2025
Version: 1.0.0
"""

from typing import List, Dict, Tuple, Optional, Any
import datetime
import os
import glob
import yaml
import pyart
import numpy as np
import gzip
import netCDF4 as nc4
from kepler_utils import (
    read_mira35_mmclx, 
    multi_mmclx2cfrad,
    cfradial_add_ncas_metadata,
    update_history_attribute,
    find_mmclxfiles,
    find_mmclx_rhi_files,
    find_mmclx_ppi_files,
    find_mmclx_vad_files,
    split_monotonic_sequence
)

def get_file_elevation(filepath: str, gzip_flag: bool = True) -> Optional[float]:
    """
    Get the mean elevation angle from an mmclx file.
    
    Args:
        filepath: Path to mmclx file
        gzip_flag: Whether file is gzip compressed
        
    Returns:
        Mean elevation angle in degrees, or None if unable to read
    """
    try:
        if gzip_flag:
            with gzip.open(filepath, 'rb') as gz:
                with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                    if 'elv' in nc.variables:
                        return float(nc.variables['elv'][:].mean())
        else:
            with nc4.Dataset(filepath, 'r') as nc:
                if 'elv' in nc.variables:
                    return float(nc.variables['elv'][:].mean())
    except Exception as e:
        print(f"Warning: Could not read elevation from {filepath}: {e}")
    return None

def apply_phase_sweeps_to_radar(radar: pyart.core.Radar, phases: List[Dict[str, Any]]) -> pyart.core.Radar:
    """
    Modify a radar object to split rays into sweeps based on detected phases.
    
    This function updates the radar object's sweep metadata to create separate
    sweeps for each detected phase. Phase types are mapped to CF-Radial 1.4
    compliant sweep_mode values:
    
    - upward_rhi / downward_rhi → 'rhi' (elevation scanning, constant azimuth)
    - dwell → 'pointing' (ALL rays elevation ≥88°, constant within 0.2°)
    - vertical_pointing → 'vertical_pointing' (ALL rays elevation ≥88°, mean >89.5°, constant within 0.2°)
    - turning → 'manual_ppi' (azimuth changing)
    - tracking / other → 'manual_rhi' (general fallback)
    
    Original phase types are preserved in metadata['phase_sequence'] for visualization.
    
    Also adds CF-Radial standard antenna_transition coordinate variable (1D, indexed by time)
    to flag individual rays with rapid angular changes (antenna repositioning between phases).
    
    Args:
        radar: PyART Radar object to modify
        phases: List of phase dictionaries from detect_man_sweep_phases()
        
    Returns:
        Modified radar object with phase-based sweeps and CF-compliant sweep_mode
    """
    n_sweeps = len(phases)
    
    # Map phase types to CF-Radial compliant sweep_mode values
    phase_to_sweep_mode = {
        'upward_rhi': 'manual_rhi',       # Manual tracking - elevation scanning upward, azimuth relatively constant
        'downward_rhi': 'manual_rhi',     # Manual tracking - elevation scanning downward, azimuth relatively constant
        'dwell': 'pointing',              # Near-vertical pointing, minimal movement
        'vertical_pointing': 'vertical_pointing',  # Near-vertical (>89.5deg), minimal movement
        'turning': 'manual_ppi',          # Azimuth changing (aircraft turning)
        'tracking': 'manual_rhi',         # General tracking fallback
        'other': 'manual_rhi'             # Ambiguous cases
    }
    
    # Create sweep metadata arrays
    sweep_number = np.arange(n_sweeps, dtype=np.int32)
    sweep_mode = np.array([phase_to_sweep_mode.get(phase['phase'], 'manual_rhi') 
                           for phase in phases], dtype='S32')
    sweep_start_ray_index = np.array([phase['start_idx'] for phase in phases], dtype=np.int32)
    sweep_end_ray_index = np.array([phase['end_idx'] for phase in phases], dtype=np.int32)
    fixed_angle = np.array([phase['mean_el'] for phase in phases], dtype=np.float32)
    
    # Update radar sweep fields - need to update the entire dictionary structure, not just data
    # to avoid dimension mismatches when changing from 1 sweep to many sweeps
    radar.sweep_number = {
        'data': sweep_number,
        'long_name': 'sweep_number',
        'units': 'count'
    }
    radar.sweep_mode = {
        'data': sweep_mode,
        'long_name': 'sweep_mode', 
        'units': 'unitless'
    }
    radar.sweep_start_ray_index = {
        'data': sweep_start_ray_index,
        'long_name': 'index_of_first_ray_in_sweep',
        'units': 'count'
    }
    radar.sweep_end_ray_index = {
        'data': sweep_end_ray_index,
        'long_name': 'index_of_last_ray_in_sweep',
        'units': 'count'
    }
    radar.fixed_angle = {
        'data': fixed_angle,
        'long_name': 'target_fixed_angle',
        'units': 'degrees'
    }
    
    # Update target_scan_rate - needs to be sized for number of sweeps
    # Calculate a representative scan rate for each phase based on its rays
    if hasattr(radar, 'target_scan_rate') and radar.target_scan_rate is not None:
        # Create target_scan_rate array for all sweeps with a default value
        target_scan_rates = np.full(n_sweeps, 1.0, dtype=np.float32)
        radar.target_scan_rate = {
            'data': target_scan_rates,
            'long_name': 'target_scan_rate_for_sweep',
            'units': 'degrees_per_second'
        }
    
    # Update scan_rate if it exists (time-dimensional, not sweep-dimensional)
    # This stays as-is since it's indexed by time, not sweep
    
    # Remove any other sweep-dimensional attributes that might exist from the original single-sweep setup
    # These need to be removed or resized to avoid dimension mismatch errors
    attrs_to_check = ['rays_per_sweep', 'sweep_type']
    for attr_name in attrs_to_check:
        if hasattr(radar, attr_name):
            attr = getattr(radar, attr_name)
            if isinstance(attr, dict) and 'data' in attr:
                # Remove or resize if it's sweep-dimensional
                if hasattr(attr['data'], 'shape') and len(attr['data'].shape) > 0:
                    if attr['data'].shape[0] == 1:  # Was sized for 1 sweep
                        # Remove it - PyART will recreate if needed
                        delattr(radar, attr_name)
                        print(f"  Removed {attr_name} attribute (was sized for 1 sweep)")
    
    # Store original phase types as metadata for visualization (non-CF-standard but useful)
    # This preserves distinction between upward_rhi/downward_rhi for quicklooks
    phase_types = [phase['phase'] for phase in phases]
    radar.metadata['phase_sequence'] = ', '.join(phase_types)
    
    # Update nsweeps attribute
    radar.nsweeps = n_sweeps
    
    # Update sweep_number per ray (assign each ray to its sweep)
    ray_sweep_number = np.zeros(radar.nrays, dtype=np.int32)
    for i, phase in enumerate(phases):
        ray_sweep_number[phase['start_idx']:phase['end_idx']+1] = i
    
    if 'sweep_number' in radar.fields:
        radar.fields['sweep_number']['data'] = ray_sweep_number
    
    # Add CF-Radial standard antenna_transition coordinate variable (not a field)
    # Mark rays with rapid angular changes (antenna repositioning between phases)
    # Also mark rays from merged 'other' phases
    antenna_transition_flags = _compute_antenna_transition_flag(radar, phases)
    
    # Set as radar attribute (1D coordinate variable, not 2D field)
    antenna_transition_dict = {
        'data': antenna_transition_flags,
        'long_name': 'antenna_is_in_transition_status',
        'units': '1',
        'comments': 'Flag indicating antenna is in transition between sweeps (1=transition, 0=stable)',
        'flag_values': np.array([0, 1], dtype=np.int8),
        'flag_meanings': 'not_in_transition in_transition'
    }
    
    # Store as radar-level coordinate (not in fields dictionary)
    radar.antenna_transition = antenna_transition_dict
    
    return radar


def detect_man_sweep_phases(azimuth: np.ndarray, elevation: np.ndarray, 
                            time: np.ndarray, 
                            elev_rate_threshold: float = 0.005,
                            az_stable_threshold: float = 0.08,
                            dwell_elevation: float = 88.0,
                            min_phase_rays: int = 10) -> List[Dict[str, Any]]:
    """
    Detect sweep phases in MAN (manual tracking) scans based on azimuth/elevation patterns.
    
    MAN scans are single files per aircraft track containing multiple phases:
    - Upward RHI: elevation increasing, azimuth relatively constant
      * Includes manual tracking (variable rate as aircraft approaches)
    - Dwell/Pointing: ALL rays with elevation ≥88°, elevation constant within 0.2° (std and range)
      * Classified as 'dwell' with sweep_mode='pointing'
      * If elevation >89.5°, classified as 'vertical_pointing'
      * Strictly validated: min elevation >= 88°, std <= 0.2°, range <= 0.2°
    - Downward RHI: elevation decreasing, azimuth relatively constant  
    - Turning: azimuth changing significantly (aircraft turning maneuvers)
    - Other: ambiguous cases (rare, expected to be minimal)
    
    Manual tracking RHIs have highly variable elevation rates (very slow when aircraft
    is far, accelerating as it approaches), so detection uses minimal directional 
    threshold (0.005 deg/s ≈ 0.3 deg/min) to catch the entire tracking sequence.
    
    This function analyzes the time series of azimuth and elevation to identify
    transition points between these phases. Small "other" segments are automatically
    merged into adjacent non-dwell phases to prevent over-fragmentation. Gaps are
    not merged into dwells unless the elevation matches within 0.2 degrees. Dwells
    are re-validated after merging to ensure strict criteria are maintained.
    
    Args:
        azimuth: Array of azimuth angles (degrees)
        elevation: Array of elevation angles (degrees)
        time: Array of time values (for computing rates of change)
        elev_rate_threshold: Minimum elevation change rate to detect motion (deg/s, default: 0.005)
        az_stable_threshold: Azimuth rate threshold for "stable" (deg/s, default: 0.08)
        dwell_elevation: Minimum elevation for vertical dwell (degrees, default: 88.0)
        min_phase_rays: Minimum rays to qualify as a phase (default: 10)
        
    Returns:
        List of dictionaries, each containing:
            - 'phase': str ('upward_rhi', 'dwell', 'vertical_pointing', 'downward_rhi', 'turning', 'other')
            - 'start_idx': int (starting ray index)
            - 'end_idx': int (ending ray index)
            - 'n_rays': int (number of rays in phase)
            - 'mean_az': float (mean azimuth for this phase)
            - 'mean_el': float (mean elevation for this phase)
            - 'std_el': float (standard deviation of elevation within phase)
            
    Example:
        >>> phases = detect_man_sweep_phases(az, el, time)
        >>> for phase in phases:
        ...     print(f"{phase['phase']}: {phase['n_rays']} rays at az={phase['mean_az']:.1f}°")
    """
    n_rays = len(elevation)
    
    if n_rays < 2:
        # Not enough rays to detect phases
        return [{
            'phase': 'tracking',
            'start_idx': 0,
            'end_idx': n_rays - 1,
            'mean_az': float(np.mean(azimuth)),
            'mean_el': float(np.mean(elevation)),
            'transition_ray_indices': []
        }]
    
    # Compute time differences (handle potential non-uniform sampling)
    dt = np.diff(time)
    dt = np.where(dt > 0, dt, 1.0)  # Avoid division by zero
    
    # Compute angular rates
    d_elv = np.diff(elevation) / dt  # deg/s
    d_az = np.diff(azimuth) / dt     # deg/s
    
    # Handle azimuth wrap-around (e.g., 359° -> 1°)
    d_az = np.where(d_az > 180, d_az - 360, d_az)
    d_az = np.where(d_az < -180, d_az + 360, d_az)
    
    # Classify each ray (using padded arrays to match original length)
    # Pad with zeros at end to match ray count
    d_elv_padded = np.concatenate([d_elv, [0]])
    d_az_padded = np.concatenate([d_az, [0]])
    
    # Define phase classifications with clear priorities
    # For manual tracking RHIs, elevation rate varies (slow when far, fast when close)
    # So we use directional change rather than fixed rate threshold
    
    # 1. Turning: azimuth changing significantly (regardless of elevation)
    #    This catches all aircraft turns, including those with varying elevation
    #    BUT: exclude high-elevation near-vertical pointing where azimuth changes don't matter
    is_turning = (np.abs(d_az_padded) > az_stable_threshold) & (elevation < dwell_elevation - 5.0)
    
    # 2. Upward RHI: elevation consistently increasing, azimuth stable
    #    Use very low threshold (0.005) to catch very slow tracking when aircraft is far
    #    This is ~0.3 deg/minute, essentially any detectable upward movement
    is_upward = (d_elv_padded > 0.005) & (np.abs(d_az_padded) < az_stable_threshold)
    
    # 3. Downward RHI: elevation consistently decreasing, azimuth stable
    is_downward = (d_elv_padded < -0.005) & (np.abs(d_az_padded) < az_stable_threshold)
    
    # 4. Dwell: minimal elevation movement at high elevation (will be refined to groups with constant elevation)
    #    Initial detection based on rate, will be post-processed to check elevation variance
    #    Azimuth can change during vertical pointing (aircraft circling overhead)
    is_dwell_candidate = (np.abs(d_elv_padded) < 0.5) & (elevation >= dwell_elevation - 2.0)
    
    # Create phase labels (priority: turning > upward > downward > dwell > other)
    phase_labels = np.full(n_rays, 'other', dtype='U20')
    phase_labels[is_dwell_candidate] = 'dwell'
    phase_labels[is_downward] = 'downward_rhi'
    phase_labels[is_upward] = 'upward_rhi'
    phase_labels[is_turning] = 'turning'
    
    # Find phase transitions
    phase_changes = np.where(phase_labels[1:] != phase_labels[:-1])[0] + 1
    phase_boundaries = np.concatenate([[0], phase_changes, [n_rays]])
    
    # Build initial phase list
    phases = []
    for i in range(len(phase_boundaries) - 1):
        start_idx = phase_boundaries[i]
        end_idx = phase_boundaries[i + 1] - 1
        phase_type = phase_labels[start_idx]
        
        phase_elevation = elevation[start_idx:end_idx+1]
        mean_el = float(np.mean(phase_elevation))
        std_el = float(np.std(phase_elevation))
        min_el = float(np.min(phase_elevation))
        max_el = float(np.max(phase_elevation))
        
        # Refine dwell classification: check elevation variance AND that all rays are near-vertical
        if phase_type == 'dwell':
            # Dwells must have:
            # 1. Reasonably constant elevation (std_el <= 1.0 deg allows for some wobble)
            # 2. Most rays at high elevation (min_el >= 85.0, allowing for edge transitions)
            # 3. Elevation range within reasonable bounds (< 3 deg)
            if std_el > 1.0 or min_el < 85.0 or (max_el - min_el) > 3.0:
                # Not a true dwell - reclassify based on context
                phase_type = 'other'
            elif mean_el > 89.5:
                # High elevation dwell -> vertical_pointing
                phase_type = 'vertical_pointing'
        
        phases.append({
            'phase': phase_type,
            'start_idx': int(start_idx),
            'end_idx': int(end_idx),
            'n_rays': int(end_idx - start_idx + 1),
            'mean_az': float(np.mean(azimuth[start_idx:end_idx+1])),
            'mean_el': mean_el,
            'std_el': std_el,
            'transition_ray_indices': []  # Track rays from merged 'other' phases
        })
    
    # Post-processing: merge ALL "other" segments into adjacent phases
    # "Other" phases are ambiguous and should not exist as standalone sweeps
    # Mark all their rays as antenna transitions
    merged_phases = []
    i = 0
    while i < len(phases):
        current = phases[i]
        
        # ALWAYS merge "other" phases into adjacent phases (regardless of size)
        # These are ambiguous classifications that shouldn't be standalone sweeps
        if current['phase'] == 'other':
            # Only mark as antenna_transition if segment is >= 30 rays (significant repositioning)
            # Smaller segments are likely just brief scanning irregularities, not true transitions
            other_ray_indices = list(range(current['start_idx'], current['end_idx'] + 1)) if current['n_rays'] >= 30 else []
            
            merged = False
            # Try to merge with previous phase first
            if merged_phases:
                prev_phase = merged_phases[-1]
                # Can always merge into non-dwell phases
                if prev_phase['phase'] not in ['dwell', 'vertical_pointing']:
                    merged_phases[-1]['end_idx'] = current['end_idx']
                    merged_phases[-1]['n_rays'] += current['n_rays']
                    merged_phases[-1]['transition_ray_indices'].extend(other_ray_indices)
                    merged = True
                # Can merge into dwell if elevation is compatible
                elif abs(current['mean_el'] - prev_phase['mean_el']) <= 0.2:
                    merged_phases[-1]['end_idx'] = current['end_idx']
                    merged_phases[-1]['n_rays'] += current['n_rays']
                    merged_phases[-1]['transition_ray_indices'].extend(other_ray_indices)
                    merged = True
            
            # If couldn't merge with previous, must merge with next
            if not merged and i < len(phases) - 1:
                next_phase = phases[i + 1]
                # Extend next phase backward to absorb this "other" phase
                phases[i + 1]['start_idx'] = current['start_idx']
                phases[i + 1]['n_rays'] += current['n_rays']
                phases[i + 1]['transition_ray_indices'].extend(other_ray_indices)
                merged = True
            
            # If still not merged (isolated other at end), convert to transition
            if not merged:
                # Convert to the most likely phase type and mark all as transitions
                if merged_phases:
                    # Adopt the previous phase type
                    current['phase'] = merged_phases[-1]['phase']
                else:
                    # No context, make it turning (safest assumption)
                    current['phase'] = 'turning'
                current['transition_ray_indices'] = other_ray_indices
                merged_phases.append(current)
            
            i += 1
            continue
        
        # Also merge very short turning phases into adjacent phases
        if current['phase'] == 'turning' and current['n_rays'] < min_phase_rays * 2:
            # Only mark as transitions if >= 30 rays (significant antenna movement)
            # Smaller segments are likely part of continuous tracking, not true transitions
            turning_ray_indices = list(range(current['start_idx'], current['end_idx'] + 1)) if current['n_rays'] >= 30 else []
            
            merged = False
            if merged_phases:
                prev_phase = merged_phases[-1]
                if prev_phase['phase'] not in ['dwell', 'vertical_pointing']:
                    merged_phases[-1]['end_idx'] = current['end_idx']
                    merged_phases[-1]['n_rays'] += current['n_rays']
                    merged_phases[-1]['transition_ray_indices'].extend(turning_ray_indices)
                    merged = True
            
            if not merged and i < len(phases) - 1:
                next_phase = phases[i + 1]
                if next_phase['phase'] not in ['dwell', 'vertical_pointing']:
                    phases[i + 1]['start_idx'] = current['start_idx']
                    phases[i + 1]['n_rays'] += current['n_rays']
                    phases[i + 1]['transition_ray_indices'].extend(turning_ray_indices)
                    merged = True
            
            if not merged:
                # Keep as turning but mark all rays as transitions
                current['transition_ray_indices'] = turning_ray_indices
                merged_phases.append(current)
            
            i += 1
            continue
        
        # Handle very short phases by merging into adjacent phases
        if current['n_rays'] < min_phase_rays and len(phases) > 1:
            # Try to merge with adjacent phase, but don't merge into dwells
            # unless elevation matches
            if merged_phases:
                prev_phase = merged_phases[-1]
                # Don't merge into dwell/vertical_pointing unless elevation matches
                if prev_phase['phase'] in ['dwell', 'vertical_pointing']:
                    # Check elevation compatibility
                    if abs(current['mean_el'] - prev_phase['mean_el']) <= 0.2:
                        # Compatible - merge into dwell
                        merged_phases[-1]['end_idx'] = current['end_idx']
                        merged_phases[-1]['n_rays'] += current['n_rays']
                        # Preserve transition ray indices if merging
                        if 'transition_ray_indices' in current:
                            merged_phases[-1]['transition_ray_indices'].extend(current['transition_ray_indices'])
                    else:
                        # Not compatible - keep separate or merge forward
                        merged_phases.append(current)
                else:
                    # Not a dwell - safe to merge into previous phase
                    merged_phases[-1]['end_idx'] = current['end_idx']
                    merged_phases[-1]['n_rays'] += current['n_rays']
                    # Preserve transition ray indices if merging
                    if 'transition_ray_indices' in current:
                        merged_phases[-1]['transition_ray_indices'].extend(current['transition_ray_indices'])
            elif i < len(phases) - 1:
                # No previous phase, keep short first phase
                merged_phases.append(current)
            else:
                # Only phase - keep it
                merged_phases.append(current)
        else:
            # Normal-sized phase or last resort
            merged_phases.append(current)
        
        i += 1
    
    # If no phases survived filtering, return single tracking phase
    if not merged_phases:
        merged_phases = [{
            'phase': 'tracking',
            'start_idx': 0,
            'end_idx': n_rays - 1,
            'n_rays': n_rays,
            'mean_az': float(np.mean(azimuth)),
            'mean_el': float(np.mean(elevation)),
            'transition_ray_indices': []
        }]
    
    # Post-processing: merge consecutive upward_rhi or downward_rhi segments
    # This joins fragmented manual RHIs that were interrupted by brief pauses or noise
    # Also extends RHIs backward to include their low-elevation starts
    final_phases = []
    skip_indices = set()  # Track phases that have been merged into others
    i = 0
    while i < len(merged_phases):
        # Skip if this phase was already merged into a previous one
        if i in skip_indices:
            i += 1
            continue
            
        current = merged_phases[i]
        
        # Try to merge consecutive RHI segments of the same direction
        if current['phase'] in ['upward_rhi', 'downward_rhi']:
            target_phase = current['phase']
            merge_start_idx = current['start_idx']
            merge_end_idx = current['end_idx']
            merge_n_rays = current['n_rays']
            all_transition_rays = current.get('transition_ray_indices', [])
            
            # FIRST: Try to extend backward to capture any preceding "other"/"turning" rays
            # that should be part of this RHI (e.g., low-elevation start with slow rate)
            if i > 0 and i-1 not in skip_indices:
                prev_phase = merged_phases[i-1]
                if prev_phase['phase'] in ['other', 'turning']:
                    # Check if this is a small segment that should be part of the RHI start
                    if prev_phase['n_rays'] < min_phase_rays * 2:  # Generous threshold
                        # For upward RHI, check if previous segment is at lower elevation
                        if target_phase == 'upward_rhi' and prev_phase['mean_el'] <= current['mean_el']:
                            # Extend start backward to include this segment
                            merge_start_idx = prev_phase['start_idx']
                            merge_n_rays = merge_end_idx - merge_start_idx + 1
                            skip_indices.add(i-1)  # Mark for skipping
                            all_transition_rays.extend(prev_phase.get('transition_ray_indices', []))
                        # For downward RHI, check if previous segment is at higher elevation  
                        elif target_phase == 'downward_rhi' and prev_phase['mean_el'] >= current['mean_el']:
                            merge_start_idx = prev_phase['start_idx']
                            merge_n_rays = merge_end_idx - merge_start_idx + 1
                            skip_indices.add(i-1)
                            all_transition_rays.extend(prev_phase.get('transition_ray_indices', []))
            
            # Look ahead for more segments of the same type
            j = i + 1
            while j < len(merged_phases):
                next_phase = merged_phases[j]
                
                # If we find another segment of the same RHI type
                if next_phase['phase'] == target_phase:
                    # Check if azimuths are similar (within 3 degrees)
                    # This indicates a continuous RHI at the same pointing direction
                    az_diff = abs(next_phase['mean_az'] - azimuth[merge_start_idx:merge_end_idx+1].mean())
                    if az_diff > 180:
                        az_diff = 360 - az_diff  # Handle wrap-around
                    
                    if az_diff <= 3.0:
                        # Similar azimuth - merge even if there's a larger gap
                        # Check the gap between segments
                        gap_start = merge_end_idx + 1
                        gap_end = next_phase['start_idx'] - 1
                        gap_size = gap_end - gap_start + 1
                        
                        # Merge if gap is small OR if gap is reasonable for continuous scan (< 50 rays)
                        if gap_size < min_phase_rays * 5:  # Allow gaps up to 50 rays
                            # Extend the merge to include the gap AND the next RHI segment
                            merge_end_idx = next_phase['end_idx']
                            merge_n_rays = merge_end_idx - merge_start_idx + 1
                            skip_indices.add(j)
                            all_transition_rays.extend(next_phase.get('transition_ray_indices', []))
                            j += 1
                            continue
                        else:
                            # Gap too large, stop merging
                            break
                    else:
                        # Different azimuth - likely a different RHI, stop merging
                        break
                        
                # If gap segment is "other" or turning (but NOT dwell/vertical_pointing)
                # Be more aggressive: merge larger gaps if they're between same-type RHIs at similar azimuth
                elif next_phase['phase'] in ['other', 'turning']:
                    # Check if there's another RHI segment of the same type after this gap
                    if j + 1 < len(merged_phases) and merged_phases[j + 1]['phase'] == target_phase:
                        following_phase = merged_phases[j + 1]
                        # Check if following phase has similar azimuth
                        az_diff = abs(following_phase['mean_az'] - azimuth[merge_start_idx:merge_end_idx+1].mean())
                        if az_diff > 180:
                            az_diff = 360 - az_diff
                        
                        # If following phase is at similar azimuth and gap isn't huge, merge through
                        if az_diff <= 3.0 and next_phase['n_rays'] < min_phase_rays * 5:  # Up to 50 rays
                            # Merge the gap phase into the RHI (keep as normal data, not transitions)
                            merge_end_idx = next_phase['end_idx']
                            merge_n_rays = merge_end_idx - merge_start_idx + 1
                            skip_indices.add(j)
                            all_transition_rays.extend(next_phase.get('transition_ray_indices', []))
                            j += 1
                            continue
                    # Also merge small gaps unconditionally (original behavior)
                    elif next_phase['n_rays'] < min_phase_rays:
                        merge_end_idx = next_phase['end_idx']
                        merge_n_rays = merge_end_idx - merge_start_idx + 1
                        skip_indices.add(j)
                        all_transition_rays.extend(next_phase.get('transition_ray_indices', []))
                        j += 1
                        continue
                    
                    # Gap too large or no following RHI, stop merging
                    break
                else:
                    # Different phase type and not a small gap, stop merging
                    break
            
            final_phases.append({
                'phase': target_phase,
                'start_idx': merge_start_idx,
                'end_idx': merge_end_idx,
                'n_rays': merge_n_rays,
                'mean_az': float(np.mean(azimuth[merge_start_idx:merge_end_idx+1])),
                'mean_el': float(np.mean(elevation[merge_start_idx:merge_end_idx+1])),
                'transition_ray_indices': all_transition_rays
            })
            
            # Move to next unprocessed phase
            i += 1
        elif current['phase'] in ['dwell', 'vertical_pointing']:
            # Merge consecutive dwell/vertical_pointing segments
            # These can be fragmented if there are brief gaps or small elevation changes
            target_phase = current['phase']
            merge_start_idx = current['start_idx']
            merge_end_idx = current['end_idx']
            merge_n_rays = current['n_rays']
            all_transition_rays = current.get('transition_ray_indices', [])
            
            # Look ahead for more dwell/vertical_pointing segments
            j = i + 1
            while j < len(merged_phases):
                next_phase = merged_phases[j]
                
                # If we find another dwell or vertical_pointing segment
                if next_phase['phase'] in ['dwell', 'vertical_pointing']:
                    # Check if elevations are compatible (within 3 degrees)
                    if abs(next_phase['mean_el'] - current['mean_el']) <= 3.0:
                        # Check the gap between segments
                        gap_start = merge_end_idx + 1
                        gap_end = next_phase['start_idx'] - 1
                        gap_size = gap_end - gap_start + 1
                        
                        # If gap is small, merge across it
                        if gap_size < min_phase_rays:
                            merge_end_idx = next_phase['end_idx']
                            merge_n_rays = merge_end_idx - merge_start_idx + 1
                            skip_indices.add(j)
                            all_transition_rays.extend(next_phase.get('transition_ray_indices', []))
                            j += 1
                            continue
                        else:
                            # Gap too large, stop merging
                            break
                    else:
                        # Elevation mismatch, stop merging
                        break
                        
                # If gap segment is small "other" or turning at similar elevation
                elif next_phase['n_rays'] < min_phase_rays and next_phase['phase'] in ['other', 'turning']:
                    # Check if this gap is at similar elevation (small drop-out)
                    if abs(next_phase['mean_el'] - current['mean_el']) <= 3.0:
                        # Merge it into the dwell to bridge the gap
                        merge_end_idx = next_phase['end_idx']
                        merge_n_rays = merge_end_idx - merge_start_idx + 1
                        skip_indices.add(j)
                        all_transition_rays.extend(next_phase.get('transition_ray_indices', []))
                        j += 1
                        continue
                    else:
                        # Different elevation, stop merging
                        break
                else:
                    # Different phase type, stop merging
                    break
            
            # Determine final phase type based on mean elevation of merged segment
            merged_elev = elevation[merge_start_idx:merge_end_idx+1]
            mean_merged_el = float(np.mean(merged_elev))
            final_phase_type = 'vertical_pointing' if mean_merged_el > 89.5 else 'dwell'
            
            final_phases.append({
                'phase': final_phase_type,
                'start_idx': merge_start_idx,
                'end_idx': merge_end_idx,
                'n_rays': merge_n_rays,
                'mean_az': float(np.mean(azimuth[merge_start_idx:merge_end_idx+1])),
                'mean_el': mean_merged_el,
                'transition_ray_indices': all_transition_rays
            })
            
            # Move to next unprocessed phase
            i += 1
        else:
            # Not an RHI phase, keep as-is (preserve transition_ray_indices)
            if 'transition_ray_indices' not in current:
                current['transition_ray_indices'] = []
            final_phases.append(current)
            i += 1
    
    # Recalculate mean angles for final phases and re-validate dwells
    for phase in final_phases:
        start = phase['start_idx']
        end = phase['end_idx']
        phase['mean_az'] = float(np.mean(azimuth[start:end+1]))
        phase['mean_el'] = float(np.mean(elevation[start:end+1]))
        
        # Re-validate dwells and vertical_pointing after merging
        if phase['phase'] in ['dwell', 'vertical_pointing']:
            phase_elevation = elevation[start:end+1]
            std_el = float(np.std(phase_elevation))
            min_el = float(np.min(phase_elevation))
            max_el = float(np.max(phase_elevation))
            
            # Use relaxed criteria to match initial detection
            if std_el > 1.0 or min_el < 85.0 or (max_el - min_el) > 3.0:
                # Invalidate - reclassify as turning (NOT "other" to avoid creating new short phases)
                # Mark all rays as transitions since this failed dwell validation
                phase['phase'] = 'turning'
                if 'transition_ray_indices' not in phase:
                    phase['transition_ray_indices'] = []
                phase['transition_ray_indices'].extend(range(start, end + 1))
            elif phase['mean_el'] > 89.5 and phase['phase'] == 'dwell':
                # Upgrade to vertical_pointing if high enough
                phase['phase'] = 'vertical_pointing'
            elif phase['mean_el'] <= 89.5 and phase['phase'] == 'vertical_pointing':
                # Downgrade to dwell if not high enough
                phase['phase'] = 'dwell'
    
    # Validate that all rays are covered with no gaps
    covered_rays = set()
    for phase in final_phases:
        phase_rays = set(range(phase['start_idx'], phase['end_idx'] + 1))
        if covered_rays & phase_rays:
            print("WARNING: Overlapping ray indices detected in phases!")
        covered_rays.update(phase_rays)
    
    all_rays = set(range(n_rays))
    if covered_rays != all_rays:
        missing_rays = all_rays - covered_rays
        print(f"WARNING: {len(missing_rays)} rays not assigned to any phase!")
        print(f"Missing ray indices: {sorted(list(missing_rays))[:20]}...")  # Show first 20
        
        # Fix by extending adjacent phases to cover gaps
        if final_phases:
            final_phases = _fill_phase_gaps(final_phases, n_rays)
            
            # Re-validate dwells after gap filling (in case any were extended)
            for phase in final_phases:
                if phase['phase'] in ['dwell', 'vertical_pointing']:
                    start = phase['start_idx']
                    end = phase['end_idx']
                    phase_elevation = elevation[start:end+1]
                    std_el = float(np.std(phase_elevation))
                    min_el = float(np.min(phase_elevation))
                    max_el = float(np.max(phase_elevation))
                    
                    # Use relaxed criteria
                    if std_el > 1.0 or min_el < 85.0 or (max_el - min_el) > 3.0:
                        # Reclassify as turning (NOT "other") and mark as transitions
                        phase['phase'] = 'turning'
                        if 'transition_ray_indices' not in phase:
                            phase['transition_ray_indices'] = []
                        phase['transition_ray_indices'].extend(range(start, end + 1))
    
    # Final safety pass: merge any remaining "other" phases that slipped through
    # This ensures no standalone "other" sweeps exist
    i = 0
    while i < len(final_phases):
        if final_phases[i]['phase'] == 'other':
            other_rays = list(range(final_phases[i]['start_idx'], final_phases[i]['end_idx'] + 1))
            
            # Try to merge with previous phase
            if i > 0:
                final_phases[i-1]['end_idx'] = final_phases[i]['end_idx']
                final_phases[i-1]['n_rays'] = final_phases[i-1]['end_idx'] - final_phases[i-1]['start_idx'] + 1
                if 'transition_ray_indices' not in final_phases[i-1]:
                    final_phases[i-1]['transition_ray_indices'] = []
                final_phases[i-1]['transition_ray_indices'].extend(other_rays)
                final_phases.pop(i)
                continue
            # Try to merge with next phase
            elif i < len(final_phases) - 1:
                final_phases[i+1]['start_idx'] = final_phases[i]['start_idx']
                final_phases[i+1]['n_rays'] = final_phases[i+1]['end_idx'] - final_phases[i+1]['start_idx'] + 1
                if 'transition_ray_indices' not in final_phases[i+1]:
                    final_phases[i+1]['transition_ray_indices'] = []
                final_phases[i+1]['transition_ray_indices'].extend(other_rays)
                final_phases.pop(i)
                continue
            else:
                # Isolated "other" phase - convert to turning and mark all as transitions
                final_phases[i]['phase'] = 'turning'
                if 'transition_ray_indices' not in final_phases[i]:
                    final_phases[i]['transition_ray_indices'] = []
                final_phases[i]['transition_ray_indices'].extend(other_rays)
        i += 1
    
    # Final enforcement: split RHI phases that don't have monotonic elevation changes
    # Each RHI phase must have strictly increasing (upward) or decreasing (downward) elevation
    monotonic_phases = []
    for phase in final_phases:
        if phase['phase'] in ['upward_rhi', 'downward_rhi']:
            # Check for monotonicity and split if needed
            split_phases = _split_non_monotonic_rhi(phase, elevation, azimuth, time)
            monotonic_phases.extend(split_phases)
        else:
            # Not an RHI - keep as-is
            monotonic_phases.append(phase)
    
    return monotonic_phases


def _split_non_monotonic_rhi(phase: Dict[str, Any], 
                              elevation: np.ndarray,
                              azimuth: np.ndarray,
                              time: np.ndarray) -> List[Dict[str, Any]]:
    """
    Split RHI phase at points where elevation direction reverses.
    
    RHI phases must have monotonic elevation changes. If a phase contains
    both upward and downward elevation movements, split it at the turning points.
    
    Args:
        phase: Phase dictionary with start_idx, end_idx, and phase type
        elevation: Full elevation array for all rays
        azimuth: Full azimuth array for all rays
        time: Full time array for all rays
        
    Returns:
        List of phase dictionaries, each with monotonic elevation
    """
    start_idx = phase['start_idx']
    end_idx = phase['end_idx']
    expected_direction = phase['phase']  # 'upward_rhi' or 'downward_rhi'
    
    # Extract elevation for this phase
    phase_elev = elevation[start_idx:end_idx+1]
    
    # Check if elevation is monotonic
    is_monotonic = True
    if expected_direction == 'upward_rhi':
        # Should be non-decreasing (allowing for flat sections)
        is_monotonic = np.all(np.diff(phase_elev) >= -0.1)  # Allow small noise
    else:  # downward_rhi
        # Should be non-increasing
        is_monotonic = np.all(np.diff(phase_elev) <= 0.1)
    
    # If monotonic, return unchanged
    if is_monotonic:
        return [phase]
    
    # Find turning points where elevation direction reverses
    elev_diff = np.diff(phase_elev)
    
    # Smooth differences to avoid splitting on noise
    # Use a simple moving average over 5 rays
    if len(elev_diff) >= 5:
        kernel = np.ones(5) / 5
        elev_diff_smooth = np.convolve(elev_diff, kernel, mode='same')
    else:
        elev_diff_smooth = elev_diff
    
    # Identify regions of upward vs downward motion (with hysteresis to avoid noise)
    is_upward = elev_diff_smooth > 0.05  # Rising
    is_downward = elev_diff_smooth < -0.05  # Falling
    
    # Create direction labels (0=unknown, 1=upward, -1=downward)
    direction = np.zeros(len(elev_diff_smooth), dtype=int)
    direction[is_upward] = 1
    direction[is_downward] = -1
    
    # Forward fill to assign unknown regions to last known direction
    for i in range(1, len(direction)):
        if direction[i] == 0 and direction[i-1] != 0:
            direction[i] = direction[i-1]
    
    # Backward fill any remaining unknowns
    for i in range(len(direction) - 2, -1, -1):
        if direction[i] == 0 and direction[i+1] != 0:
            direction[i] = direction[i+1]
    
    # Find transitions where direction changes
    direction_changes = np.where(direction[1:] != direction[:-1])[0] + 1
    
    # If no clear transitions found, return original phase
    if len(direction_changes) == 0:
        return [phase]
    
    # Create split points (indices relative to full array)
    split_boundaries = [start_idx] + [start_idx + idx for idx in direction_changes] + [end_idx + 1]
    
    # Create new phases for each monotonic segment
    split_phases = []
    for i in range(len(split_boundaries) - 1):
        seg_start = split_boundaries[i]
        seg_end = split_boundaries[i + 1] - 1
        seg_n_rays = seg_end - seg_start + 1
        
        # Keep all segments from splitting - they were part of a valid phase
        # Don't skip short segments as they might be important (e.g., start of RHI at low elevation)
        
        # Determine actual direction for this segment
        seg_elev = elevation[seg_start:seg_end+1]
        elev_trend = seg_elev[-1] - seg_elev[0]
        
        # Use lower threshold (0.2 degrees) to capture even small elevation changes
        # This ensures we keep the low-elevation start of upward RHIs
        if elev_trend > 0.2:  # Net upward movement
            seg_phase = 'upward_rhi'
        elif elev_trend < -0.2:  # Net downward movement
            seg_phase = 'downward_rhi'
        else:
            # Very small elevation change - inherit from original phase if possible
            # This handles near-horizontal segments at the beginning/end
            if abs(elev_trend) < 0.2:
                # Use the original expected direction for flat segments
                seg_phase = expected_direction
            else:
                seg_phase = 'turning'
        
        # Mark boundary rays as transitions
        transition_rays = []
        if i > 0:  # Not the first segment - mark first few rays as transition
            transition_rays.extend(range(seg_start, min(seg_start + 3, seg_end + 1)))
        if i < len(split_boundaries) - 2:  # Not the last segment - mark last few rays
            transition_rays.extend(range(max(seg_end - 2, seg_start), seg_end + 1))
        
        split_phases.append({
            'phase': seg_phase,
            'start_idx': int(seg_start),
            'end_idx': int(seg_end),
            'n_rays': int(seg_n_rays),
            'mean_az': float(np.mean(azimuth[seg_start:seg_end+1])),
            'mean_el': float(np.mean(elevation[seg_start:seg_end+1])),
            'transition_ray_indices': transition_rays
        })
    
    # If no valid segments created, return original phase
    if not split_phases:
        return [phase]
    
    return split_phases


def _compute_antenna_transition_flag(radar: pyart.core.Radar,
                                    phases: List[Dict[str, Any]] = None,
                                    transition_threshold: float = 3.0,
                                    max_transition_duration: float = 5.0) -> np.ndarray:
    """
    Compute CF-Radial standard antenna_transition flags for each ray.
    
    Returns a 1D array of flags (one per ray) to be stored as a coordinate variable,
    not as a 2D data field.
    
    Antenna transitions occur when the radar is repositioning between different tracking
    phases (e.g., from upward RHI to vertical dwell, or from dwell to downward RHI).
    These periods have:
    - Very high angular rates (>3 deg/s in elevation or azimuth)
    - Short duration (typically < 5 seconds)
    - Occur at boundaries between major phases
    - Rays from merged 'other' phases (very short ambiguous segments)
    
    Args:
        radar: PyART Radar object
        phases: Optional list of phase dictionaries with transition_ray_indices
        transition_threshold: Minimum angular rate for transition (deg/s, default 3.0)
        max_transition_duration: Maximum duration for transition (seconds, default 5.0)
        
    Returns:
        1D array of flags (0=stable, 1=transition) for each ray
    """
    n_rays = radar.nrays
    antenna_transition = np.zeros(n_rays, dtype=np.int8)
    
    if n_rays < 2:
        return antenna_transition
    
    # Extract angle and time data
    azimuth = radar.azimuth['data']
    elevation = radar.elevation['data']
    time_data = radar.time['data']
    
    # Compute time differences and angular rates
    dt = np.diff(time_data)
    dt = np.where(dt > 0, dt, 1.0)  # Avoid division by zero
    
    d_elv = np.diff(elevation) / dt  # deg/s
    d_az = np.diff(azimuth) / dt     # deg/s
    
    # Handle azimuth wrap-around
    d_az = np.where(d_az > 180, d_az - 360, d_az)
    d_az = np.where(d_az < -180, d_az + 360, d_az)
    
    # Identify rays with very high angular rates
    high_rate = (np.abs(d_elv) > transition_threshold) | (np.abs(d_az) > transition_threshold)
    
    # Mark transitions (pad to match ray count)
    # Only mark the ray where high rate occurs, not adjacent rays
    for i in range(n_rays - 1):
        if high_rate[i]:
            antenna_transition[i] = 1
    
    # Additional check: Look for very short segments between stable phases
    # that are likely antenna repositioning (not normal scanning)
    if radar.nsweeps > 1:
        for sweep_idx in range(radar.nsweeps):
            start_idx = radar.sweep_start_ray_index['data'][sweep_idx]
            end_idx = radar.sweep_end_ray_index['data'][sweep_idx]
            n_rays_sweep = end_idx - start_idx + 1
            
            # Only mark very short sweeps (<10 rays) with predominantly high angular rates
            if n_rays_sweep < 10:
                high_rate_fraction = np.sum(antenna_transition[start_idx:end_idx+1]) / n_rays_sweep
                if high_rate_fraction > 0.5:
                    # Mark entire sweep as transition
                    antenna_transition[start_idx:end_idx+1] = 1
    
    # Mark rays that came from merged 'other' phases
    if phases is not None:
        for phase in phases:
            if 'transition_ray_indices' in phase and phase['transition_ray_indices']:
                for ray_idx in phase['transition_ray_indices']:
                    if 0 <= ray_idx < n_rays:
                        antenna_transition[ray_idx] = 1
    
    return antenna_transition


def _fill_phase_gaps(phases: List[Dict[str, Any]], n_rays: int) -> List[Dict[str, Any]]:
    """
    Ensure phases cover all ray indices by extending phases to fill gaps.
    
    Gaps are assigned to adjacent non-dwell phases. Gaps should only be assigned
    to dwells/vertical_pointing if the gap elevation matches the dwell elevation
    within 0.2 deg.
    
    Args:
        phases: List of phase dictionaries
        n_rays: Total number of rays
        
    Returns:
        Modified phases list with no gaps
    """
    if not phases:
        return phases
    
    # Sort phases by start_idx to process in order
    phases = sorted(phases, key=lambda p: p['start_idx'])
    
    # Fill gap at start (if any)
    if phases[0]['start_idx'] > 0:
        phases[0]['n_rays'] += phases[0]['start_idx']
        phases[0]['start_idx'] = 0
    
    # Fill gaps between phases
    for i in range(len(phases) - 1):
        current_end = phases[i]['end_idx']
        next_start = phases[i + 1]['start_idx']
        
        if next_start > current_end + 1:
            # Gap exists - decide which phase to extend
            current_phase = phases[i]['phase']
            next_phase = phases[i + 1]['phase']
            
            # Don't extend dwells/vertical_pointing into gaps unless explicitly appropriate
            if current_phase in ['dwell', 'vertical_pointing']:
                # Extend next phase backward instead
                gap_size = next_start - current_end - 1
                phases[i + 1]['start_idx'] = current_end + 1
                phases[i + 1]['n_rays'] += gap_size
            elif next_phase in ['dwell', 'vertical_pointing']:
                # Extend current phase forward
                gap_size = next_start - current_end - 1
                phases[i]['end_idx'] = next_start - 1
                phases[i]['n_rays'] += gap_size
            else:
                # Neither is a dwell - extend current phase forward
                gap_size = next_start - current_end - 1
                phases[i]['end_idx'] = next_start - 1
                phases[i]['n_rays'] += gap_size
    
    # Fill gap at end (if any)
    if phases[-1]['end_idx'] < n_rays - 1:
        gap_size = n_rays - 1 - phases[-1]['end_idx']
        phases[-1]['end_idx'] = n_rays - 1
        phases[-1]['n_rays'] += gap_size
    
    return phases

def process_kepler_woest_day_step1(
    datestr: str, 
    inpath: str, 
    outpath: str, 
    yaml_project_file: str, 
    yaml_instrument_file: str, 
    azimuth_offset: float = -6.85, 
    revised_northangle: float = 302.15, 
    gzip_flag: bool = False, 
    data_version: str = "1.0.0",
    tracking_tag: str = 'AMOF_20220301000000',
    **kwargs
) -> None:
    """
    Process WOEST campaign data for a single day - Step 1.
    
    This function processes scans from the WOEST campaign, converting
    mmclx files to CF-Radial format. Creates multi-sweep files by grouping:
    - HSRHI: groups of 6 files (one complete scan cycle) in 30-min windows
    - BLPPI: all files in 30-min windows
    - VPT: all files for the day
    - VAD: all files for the day
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing mmclx files (with subdirs: hsrhi, blppi, vad, vpt)
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        azimuth_offset: Azimuth offset correction (default: -6.85)
        revised_northangle: North angle correction (default: 302.15)
        gzip_flag: Whether input files are gzip compressed
        data_version: Data version string
        tracking_tag: AMOF tracking tag
    """
    print(f"Processing WOEST day: {datestr}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Define subdirectory paths for WOEST data organization
    # WOEST structure: /mom/YYYYMMDD/scantype/*.mmclx.gz
    date_inpath = os.path.join(inpath, datestr)
    hsrhipath = os.path.join(date_inpath, 'hsrhi')
    blppipath = os.path.join(date_inpath, 'blppi')
    vadpath = os.path.join(date_inpath, 'vad')
    vptpath = os.path.join(date_inpath, 'vpt')
    
    # Parse date
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d')
    end_date = start_date + datetime.timedelta(days=1)
    
    # Process HSRHI and BLPPI in 30-minute windows
    current_date = start_date
    while current_date <= end_date:
        next_halfhour = current_date + datetime.timedelta(minutes=30)
        
        # Process HSRHI sector 1: -10° to 169°
        try:
            print(f"Searching for HSRHI sector 1 files ({current_date.strftime('%H:%M')})")
            hsrhi1_files = find_mmclx_rhi_files(
                current_date.strftime('%Y-%m-%d %H:%M:%S'),
                next_halfhour.strftime('%Y-%m-%d %H:%M:%S'),
                -10, 169, hsrhipath, gzip_flag=gzip_flag,
                azimuth_offset=azimuth_offset, revised_northangle=revised_northangle
            )
            
            if hsrhi1_files:
                print(f"Found {len(hsrhi1_files)} HSRHI sector 1 files")
                # Group in sets of 6 (one complete scan cycle)
                for i in range(0, len(hsrhi1_files), 6):
                    subset = hsrhi1_files[i:i+6]
                    if subset:  # Only process if there are files
                        _process_multi_woest_files(
                            subset, outdir, yaml_project_file, yaml_instrument_file,
                            'HSRHI', revised_northangle, gzip_flag, data_version,
                            tracking_tag, azimuth_offset
                        )
        except Exception as e:
            print(f"Error processing HSRHI sector 1: {e}")
        
        # Process HSRHI sector 2: 170° to 349°
        try:
            print(f"Searching for HSRHI sector 2 files ({current_date.strftime('%H:%M')})")
            hsrhi2_files = find_mmclx_rhi_files(
                current_date.strftime('%Y-%m-%d %H:%M:%S'),
                next_halfhour.strftime('%Y-%m-%d %H:%M:%S'),
                170, 349, hsrhipath, gzip_flag=gzip_flag,
                azimuth_offset=azimuth_offset, revised_northangle=revised_northangle
            )
            
            if hsrhi2_files:
                print(f"Found {len(hsrhi2_files)} HSRHI sector 2 files")
                # Group in sets of 6 (one complete scan cycle)
                for i in range(0, len(hsrhi2_files), 6):
                    subset = hsrhi2_files[i:i+6]
                    if subset:  # Only process if there are files
                        _process_multi_woest_files(
                            subset, outdir, yaml_project_file, yaml_instrument_file,
                            'HSRHI', revised_northangle, gzip_flag, data_version,
                            tracking_tag, azimuth_offset
                        )
        except Exception as e:
            print(f"Error processing HSRHI sector 2: {e}")
        
        # Process BLPPI files (all in window together)
        try:
            print(f"Searching for BLPPI files ({current_date.strftime('%H:%M')})")
            blppi_files = find_mmclx_ppi_files(
                current_date.strftime('%Y-%m-%d %H:%M:%S'),
                next_halfhour.strftime('%Y-%m-%d %H:%M:%S'),
                0, 80, blppipath, gzip_flag=gzip_flag
            )
            
            if blppi_files:
                print(f"Found {len(blppi_files)} BLPPI files")
                _process_multi_woest_files(
                    blppi_files, outdir, yaml_project_file, yaml_instrument_file,
                    'BLPPI', revised_northangle, gzip_flag, data_version,
                    tracking_tag, azimuth_offset
                )
        except Exception as e:
            print(f"Error processing BLPPI: {e}")
        
        current_date = next_halfhour
    
    # Process VPT files for entire day
    try:
        print("Searching for VPT files (full day)")
        vpt_files = find_mmclxfiles(
            start_date.strftime('%Y-%m-%d %H:%M:%S'),
            end_date.strftime('%Y-%m-%d %H:%M:%S'),
            'vert', vptpath, gzip_flag=gzip_flag
        )
        
        if vpt_files:
            print(f"Found {len(vpt_files)} VPT files")
            _process_multi_woest_files(
                vpt_files, outdir, yaml_project_file, yaml_instrument_file,
                'VPT', revised_northangle, gzip_flag, data_version,
                tracking_tag, azimuth_offset
            )
    except Exception as e:
        print(f"Error processing VPT: {e}")
    
    # Process VAD files for entire day
    try:
        print("Searching for VAD files (full day)")
        vad_files = find_mmclx_vad_files(
            start_date.strftime('%Y-%m-%d %H:%M:%S'),
            end_date.strftime('%Y-%m-%d %H:%M:%S'),
            80, 90, vadpath, gzip_flag=gzip_flag
        )
        
        if vad_files:
            print(f"Found {len(vad_files)} VAD files")
            _process_multi_woest_files(
                vad_files, outdir, yaml_project_file, yaml_instrument_file,
                'VAD', revised_northangle, gzip_flag, data_version,
                tracking_tag, azimuth_offset
            )
    except Exception as e:
        print(f"Error processing VAD: {e}")

def _process_multi_woest_files(
    infiles: List[str],
    outdir: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    scan_name: str,
    revised_northangle: float,
    gzip_flag: bool,
    data_version: str,
    tracking_tag: str,
    azimuth_offset: float = -6.85
) -> None:
    """
    Process multiple WOEST files into a single multi-sweep CF-Radial file.
    
    Args:
        infiles: List of input mmclx file paths
        outdir: Output directory
        yaml_project_file: Project YAML file
        yaml_instrument_file: Instrument YAML file
        scan_name: Scan type ('HSRHI', 'BLPPI', 'VPT', 'VAD')
        revised_northangle: North angle correction
        gzip_flag: Whether files are gzip compressed
        data_version: Data version string
        tracking_tag: AMOF tracking tag
        azimuth_offset: Azimuth offset correction
    """
    if not infiles:
        return
    
    print(f"Processing {len(infiles)} {scan_name} files into multi-sweep file")
    
    try:
        # Use multi_mmclx2cfrad to create multi-sweep file
        radar_ds = multi_mmclx2cfrad(
            infiles,
            outdir,
            scan_name=scan_name,
            campaign='woest',
            gzip_flag=gzip_flag,
            azimuth_offset=azimuth_offset,
            revised_northangle=revised_northangle,
            data_version=data_version,
            tracking_tag=tracking_tag,
            yaml_project_file=yaml_project_file,
            yaml_instrument_file=yaml_instrument_file,
            single_sweep=False  # Multi-sweep mode
        )
        
        if radar_ds is not None:
            print(f"Successfully created multi-sweep {scan_name} file")
        else:
            print(f"Warning: multi_mmclx2cfrad returned None for {scan_name}")
            
    except Exception as e:
        print(f"Error in _process_multi_woest_files for {scan_name}: {e}")
        raise


def _process_single_woest_file(
    infile: str,
    outdir: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    revised_northangle: float,
    gzip_flag: bool,
    data_version: str,
    tracking_tag: str = 'AMOF_20220301000000'
) -> None:
    """
    Process a single WOEST file.
    
    Args:
        infile: Input mmclx file path
        outdir: Output directory
        yaml_project_file: Project YAML file
        yaml_instrument_file: Instrument YAML file
        revised_northangle: North angle correction
        gzip_flag: Whether file is gzip compressed
        data_version: Data version string
    """
    # Read radar data
    radar_ds = read_mira35_mmclx(infile, gzip_flag=gzip_flag, revised_northangle=revised_northangle)
    
    # Generate output filename
    base_filename = os.path.basename(infile).replace('.mmclx', '').replace('.gz', '')
    outfile = os.path.join(outdir, f"{base_filename}_l1_v{data_version}.nc")
    
    # Write CF-Radial file
    pyart.io.write_cfradial(outfile, radar_ds, format='NETCDF4', time_reference=True)
    
    # Add NCAS metadata
    cfradial_add_ncas_metadata(outfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version)
    
    # Update history
    update_history_attribute(outfile, f"WOEST Step 1: Convert mmclx to CF-Radial, revised_northangle={revised_northangle}")
    
    print(f"Created: {outfile}")

def process_kepler_general_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    azimuth_offset: float = 0.0,
    revised_northangle: float = 55.7,
    gzip_flag: bool = False,
    data_version: str = "1.0.0",
    single_sweep: bool = True,
    tracking_tag: str = 'AMOF_20230401000000',
    campaign: str = 'ccrest-m',
    **kwargs
) -> None:
    """
    General processor for campaigns with standard RHI/PPI/VPT structure - Step 1.
    
    This function processes RHI, PPI, and VPT scans from campaigns including
    CCREST-M and COALESC3.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing mmclx files
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        azimuth_offset: Azimuth offset correction (default: 0.0)
        revised_northangle: North angle correction (default: 55.7)
        gzip_flag: Whether input files are gzip compressed
        data_version: Data version string
        single_sweep: Create separate files for each sweep (default: True)
        tracking_tag: AMOF tracking tag
        campaign: Campaign name (e.g., 'ccrest-m', 'coalesc3')
    """
    print(f"Processing {campaign.upper()} day: {datestr}")
    print(f"Single sweep mode: {single_sweep}")
    print(f"Data version: {data_version}")
    print(f"Tracking tag: {tracking_tag}")

    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Time range for the day
    start_time = f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:8]} 00:00:00"
    end_time = f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:8]} 23:59:59"
    
    # Process VPT (vertical pointing) files if present
    print("Processing VPT files...")
    try:
        vpt_files = find_mmclxfiles(start_time, end_time, 'vert', inpath, gzip_flag=gzip_flag)
        print(f"Found {len(vpt_files)} VPT files")
        
        if len(vpt_files) > 0:
            print("VPT files will always create multi-sweep output (ignoring single_sweep flag)")
            RadarDS_VPT = multi_mmclx2cfrad(
                vpt_files, outdir, scan_name='VPT', gzip_flag=gzip_flag,
                azimuth_offset=azimuth_offset, 
                tracking_tag=tracking_tag,
                campaign=campaign,  # Use campaign parameter instead of hardcoded value
                revised_northangle=revised_northangle,
                data_version=data_version,
                single_sweep=False,  # ALWAYS FALSE FOR VPT
                yaml_project_file=yaml_project_file,
                yaml_instrument_file=yaml_instrument_file
            )
            print(f"Processed {len(vpt_files)} VPT files into multi-sweep dataset")
    except Exception as e:
        print(f"VPT processing problem: {e}")
        import traceback
        traceback.print_exc()
    
    # Process RHI files
    print("Processing RHI files...")
    rhi_files = find_mmclxfiles(start_time, end_time, "rhi", inpath, gzip_flag=gzip_flag)
    print(f"Found {len(rhi_files)} RHI files")
    
    if len(rhi_files) > 0:
        if single_sweep:
            # Process each RHI file separately
            for f in rhi_files:
                try:
                    RadarDS_RHI = multi_mmclx2cfrad(
                        [f], outdir, scan_name='RHI', gzip_flag=gzip_flag,
                        azimuth_offset=azimuth_offset,
                        tracking_tag=tracking_tag,
                        campaign=campaign,
                        revised_northangle=revised_northangle,
                        data_version=data_version,
                        single_sweep=True,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    print(f"Processed single RHI file: {f}")
                except Exception as e:
                    print(f"Error processing RHI file {f}: {e}")
        else:
            # Process all RHI files together (multi-sweep)
            try:
                RadarDS_RHI = multi_mmclx2cfrad(
                    rhi_files, outdir, scan_name='RHI', gzip_flag=gzip_flag,
                    azimuth_offset=azimuth_offset,
                    tracking_tag=tracking_tag,
                    campaign=campaign,
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(rhi_files)} RHI files into combined dataset")
            except Exception as e:
                print(f"Error processing RHI files: {e}")
    
    # Process PPI files (but check elevation - high elevation scans are VAD)
    print("Processing PPI files...")
    ppi_files = find_mmclxfiles(start_time, end_time, "ppi", inpath, gzip_flag=gzip_flag)
    print(f"Found {len(ppi_files)} PPI files")
    
    # Separate PPI files by elevation (>10° are VAD scans)
    true_ppi_files = []
    vad_from_ppi_files = []
    
    for f in ppi_files:
        elv = get_file_elevation(f, gzip_flag=gzip_flag)
        if elv is not None:
            if elv > 10.0:
                vad_from_ppi_files.append(f)
                print(f"  {os.path.basename(f)}: elv={elv:.1f}° -> VAD")
            else:
                true_ppi_files.append(f)
                print(f"  {os.path.basename(f)}: elv={elv:.1f}° -> PPI")
        else:
            # If we can't read elevation, assume it's a true PPI
            true_ppi_files.append(f)
            print(f"  {os.path.basename(f)}: elv=unknown -> PPI (default)")
    
    print(f"Classified {len(true_ppi_files)} as true PPI, {len(vad_from_ppi_files)} as VAD")
    
    # Process true PPI files
    if len(true_ppi_files) > 0:
        if single_sweep:
            # Process each PPI file separately
            for f in true_ppi_files:
                try:
                    RadarDS_PPI = multi_mmclx2cfrad(
                        [f], outdir, scan_name='PPI', gzip_flag=gzip_flag,
                        azimuth_offset=azimuth_offset,
                        tracking_tag=tracking_tag,
                        campaign=campaign,
                        revised_northangle=revised_northangle,
                        data_version=data_version,
                        single_sweep=True,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    print(f"Processed single PPI file: {f}")
                except Exception as e:
                    print(f"Error processing PPI file {f}: {e}")
        else:
            # Process all PPI files together (multi-sweep)
            try:
                RadarDS_PPI = multi_mmclx2cfrad(
                    true_ppi_files, outdir, scan_name='PPI', gzip_flag=gzip_flag,
                    azimuth_offset=azimuth_offset,
                    tracking_tag=tracking_tag,
                    campaign=campaign,
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(true_ppi_files)} PPI files into combined dataset")
            except Exception as e:
                print(f"Error processing PPI files: {e}")
    
    # Process high-elevation PPI files as VAD (group by elevation, multi-sweep per elevation)
    if len(vad_from_ppi_files) > 0:
        # Group VAD files by elevation angle
        from collections import defaultdict
        vad_by_elevation = defaultdict(list)
        
        for f in vad_from_ppi_files:
            elv = get_file_elevation(f, gzip_flag=gzip_flag)
            if elv is not None:
                # Round to nearest degree to group similar elevations
                elv_rounded = round(elv)
                vad_by_elevation[elv_rounded].append(f)
        
        print(f"Found {len(vad_by_elevation)} distinct VAD elevation angles")
        
        # Process each elevation group separately
        for elv, files in sorted(vad_by_elevation.items()):
            print(f"Processing {len(files)} VAD scans at {elv}° elevation...")
            try:
                # Include elevation in scan name to create unique filenames
                scan_name = f'vad-{elv}deg' if len(vad_by_elevation) > 1 else 'vad'
                RadarDS_VAD = multi_mmclx2cfrad(
                    files, outdir, scan_name=scan_name, gzip_flag=gzip_flag,
                    azimuth_offset=azimuth_offset,
                    tracking_tag=tracking_tag,
                    campaign=campaign,
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,  # VAD always multi-sweep (one file per elevation per day)
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(files)} VAD scans at {elv}° into single multi-sweep file")
            except Exception as e:
                print(f"Error processing VAD files at {elv}°: {e}")
                import traceback
                traceback.print_exc()
    
    print(f"Completed {campaign.upper()} processing for {datestr}")

def process_kepler_picasso_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    azimuth_offset: float = 0.0,
    revised_northangle: float = 97.1,  # PICASSO north angle (should be loaded from YAML via process_campaign_day)
    gzip_flag: bool = False,
    data_version: str = "1.0.0",
    single_sweep: bool = True,
    tracking_tag: str = 'CFARR_0002',
    campaign: str = 'picasso',
    split_man_phases: bool = True,  # Split MAN scans into phase-based sweeps (upward/dwell/downward/turning)
    **kwargs
) -> None:
    """
    PICASSO-specific processor for RHI/PPI/VPT/MAN scans - Step 1.
    
    This function processes standard scans (RHI, PPI, VPT) plus special MAN (manual)
    scans that track aircraft with simultaneous azimuth and elevation changes.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing mmclx files
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        azimuth_offset: Azimuth offset correction (default: 0.0)
        revised_northangle: North angle correction (default: 97.1 for PICASSO)
        gzip_flag: Whether input files are gzip compressed
        data_version: Data version string
        single_sweep: Create separate files for each sweep (default: True)
        tracking_tag: CFARR tracking tag
        campaign: Campaign name (default: 'picasso')
        split_man_phases: If True, split MAN scans into separate sweeps per phase (default: True)
    """
    print(f"Processing {campaign.upper()} day: {datestr}")
    print(f"Single sweep mode: {single_sweep}")
    print(f"Data version: {data_version}")
    print(f"Tracking tag: {tracking_tag}")

    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Time range for the day
    start_time = f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:8]} 00:00:00"
    end_time = f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:8]} 23:59:59"
    
    # Process VPT (vertical pointing) files if present
    print("Processing VPT files...")
    try:
        vpt_files = find_mmclxfiles(start_time, end_time, 'vert', inpath, gzip_flag=gzip_flag)
        print(f"Found {len(vpt_files)} VPT files")
        
        if len(vpt_files) > 0:
            print("VPT files will always create multi-sweep output (ignoring single_sweep flag)")
            RadarDS_VPT = multi_mmclx2cfrad(
                vpt_files, outdir, scan_name='VPT', gzip_flag=gzip_flag,
                azimuth_offset=azimuth_offset, 
                tracking_tag=tracking_tag,
                campaign=campaign,
                revised_northangle=revised_northangle,
                data_version=data_version,
                single_sweep=False,  # ALWAYS FALSE FOR VPT
                yaml_project_file=yaml_project_file,
                yaml_instrument_file=yaml_instrument_file
            )
            print(f"Processed {len(vpt_files)} VPT files into multi-sweep dataset")
    except Exception as e:
        print(f"VPT processing problem: {e}")
        import traceback
        traceback.print_exc()
    
    # Process RHI files
    print("Processing RHI files...")
    rhi_files = find_mmclxfiles(start_time, end_time, "rhi", inpath, gzip_flag=gzip_flag)
    print(f"Found {len(rhi_files)} RHI files")
    
    if len(rhi_files) > 0:
        if single_sweep:
            # Process each RHI file separately
            for f in rhi_files:
                try:
                    RadarDS_RHI = multi_mmclx2cfrad(
                        [f], outdir, scan_name='RHI', gzip_flag=gzip_flag,
                        azimuth_offset=azimuth_offset,
                        tracking_tag=tracking_tag,
                        campaign=campaign,
                        revised_northangle=revised_northangle,
                        data_version=data_version,
                        single_sweep=True,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    print(f"Processed single RHI file: {f}")
                except Exception as e:
                    print(f"Error processing RHI file {f}: {e}")
        else:
            # Process all RHI files together (multi-sweep)
            try:
                RadarDS_RHI = multi_mmclx2cfrad(
                    rhi_files, outdir, scan_name='RHI', gzip_flag=gzip_flag,
                    azimuth_offset=azimuth_offset,
                    tracking_tag=tracking_tag,
                    campaign=campaign,
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(rhi_files)} RHI files into combined dataset")
            except Exception as e:
                print(f"Error processing RHI files: {e}")
    
    # Process PPI files (but check elevation - high elevation scans are VAD)
    print("Processing PPI files...")
    ppi_files = find_mmclxfiles(start_time, end_time, "ppi", inpath, gzip_flag=gzip_flag)
    print(f"Found {len(ppi_files)} PPI files")
    
    # Separate PPI files by elevation (>10° are VAD scans)
    true_ppi_files = []
    vad_from_ppi_files = []
    
    for f in ppi_files:
        elv = get_file_elevation(f, gzip_flag=gzip_flag)
        if elv is not None:
            if elv > 10.0:
                vad_from_ppi_files.append(f)
                print(f"  {os.path.basename(f)}: elv={elv:.1f}° -> VAD")
            else:
                true_ppi_files.append(f)
                print(f"  {os.path.basename(f)}: elv={elv:.1f}° -> PPI")
        else:
            # If we can't read elevation, assume it's a true PPI
            true_ppi_files.append(f)
            print(f"  {os.path.basename(f)}: elv=unknown -> PPI (default)")
    
    print(f"Classified {len(true_ppi_files)} as true PPI, {len(vad_from_ppi_files)} as VAD")
    
    # Process true PPI files
    if len(true_ppi_files) > 0:
        if single_sweep:
            # Process each PPI file separately
            for f in true_ppi_files:
                try:
                    RadarDS_PPI = multi_mmclx2cfrad(
                        [f], outdir, scan_name='PPI', gzip_flag=gzip_flag,
                        azimuth_offset=azimuth_offset,
                        tracking_tag=tracking_tag,
                        campaign=campaign,
                        revised_northangle=revised_northangle,
                        data_version=data_version,
                        single_sweep=True,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    print(f"Processed single PPI file: {f}")
                except Exception as e:
                    print(f"Error processing PPI file {f}: {e}")
        else:
            # Process all PPI files together (multi-sweep)
            try:
                RadarDS_PPI = multi_mmclx2cfrad(
                    true_ppi_files, outdir, scan_name='PPI', gzip_flag=gzip_flag,
                    azimuth_offset=azimuth_offset,
                    tracking_tag=tracking_tag,
                    campaign=campaign,
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(true_ppi_files)} PPI files into combined dataset")
            except Exception as e:
                print(f"Error processing PPI files: {e}")
    
    # Process high-elevation PPI files as VAD (group by elevation, multi-sweep per elevation)
    if len(vad_from_ppi_files) > 0:
        # Group VAD files by elevation angle
        from collections import defaultdict
        vad_by_elevation = defaultdict(list)
        
        for f in vad_from_ppi_files:
            elv = get_file_elevation(f, gzip_flag=gzip_flag)
            if elv is not None:
                # Round to nearest degree to group similar elevations
                elv_rounded = round(elv)
                vad_by_elevation[elv_rounded].append(f)
        
        print(f"Found {len(vad_by_elevation)} distinct VAD elevation angles")
        
        # Process each elevation group separately
        for elv, files in sorted(vad_by_elevation.items()):
            print(f"Processing {len(files)} VAD scans at {elv}° elevation...")
            try:
                # Include elevation in scan name to create unique filenames
                scan_name = f'vad-{elv}deg' if len(vad_by_elevation) > 1 else 'vad'
                RadarDS_VAD = multi_mmclx2cfrad(
                    files, outdir, scan_name=scan_name, gzip_flag=gzip_flag,
                    azimuth_offset=azimuth_offset,
                    tracking_tag=tracking_tag,
                    campaign=campaign,
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,  # VAD always multi-sweep (one file per elevation per day)
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(files)} VAD scans at {elv}° into single multi-sweep file")
            except Exception as e:
                print(f"Error processing VAD files at {elv}°: {e}")
    
    # Process MAN (manual tracking) scans - PICASSO-specific
    # These track aircraft with simultaneous azimuth and elevation changes
    # 
    # MAN scan structure:
    #   - Each file = one complete aircraft track
    #   - Typical phases within each file:
    #     * Upward scanning RHI (elevation increasing, azimuth ~constant)
    #     * Vertical dwell (elevation ≥88°, both angles stable within 0.2°)
    #     * Downward scanning RHI (elevation decreasing, azimuth ~constant)
    #     * Antenna transitions (rapid repositioning, flagged with CF antenna_transition field)
    #     * Possible turning segments (both az/el changing)
    #
    # Processing approach:
    #   - Each MAN file processed as multi-sweep netCDF (one file per aircraft track)
    #   - Optional: split into separate sweeps per phase using detect_man_sweep_phases()
    #   - Phase types mapped to CF-Radial 1.4 compliant sweep_mode values:
    #     * upward_rhi/downward_rhi → 'rhi' (elevation scanning, constant azimuth)
    #     * dwell → 'pointing' (near-vertical <89.5°, constant elevation within 0.2°)
    #     * vertical_pointing → 'vertical_pointing' (near-vertical >89.5°, constant within 0.2°)
    #     * turning → 'manual_ppi' (azimuth changing)
    #   - Original phase types preserved in metadata['phase_sequence'] for visualization
    #   - Gaps only merged into dwells if elevation matches within 0.2°
    #   - CF-Radial antenna_transition coordinate variable added to flag rapid slewing between phases
    print("Processing MAN (manual tracking) files...")
    try:
        man_files = find_mmclxfiles(start_time, end_time, "man", inpath, gzip_flag=gzip_flag)
        print(f"Found {len(man_files)} MAN files")
        
        if len(man_files) > 0:
            # Get platform location from YAML file for filename construction
            location = 'unknown'  # default
            try:
                with open(yaml_project_file, 'r') as file:
                    projects = yaml.safe_load(file)
                
                # Find the project with the matching tracking_tag
                project = None
                for p in projects:
                    if tracking_tag in p:
                        project = p[tracking_tag]
                        break
                
                if project and 'ncas_instruments' in project:
                    # Find the radar instrument
                    for instrument in project['ncas_instruments']:
                        if 'ncas-mobile-ka-band-radar-1' in instrument:
                            radar_info = instrument['ncas-mobile-ka-band-radar-1']
                            if 'platform' in radar_info and 'location' in radar_info['platform']:
                                location = radar_info['platform']['location'].lower()
                                print(f"Using location from YAML: {location}")
                                break
                
                if location == 'unknown':
                    print(f"Warning: Could not find platform location in {yaml_project_file}, using campaign name")
                    location = campaign  # fallback to campaign name
            
            except Exception as e:
                print(f"Warning: Error reading location from {yaml_project_file}: {e}")
                location = campaign  # fallback to campaign name
            
            # Process each MAN file as one multi-sweep netCDF (one file per aircraft track)
            if split_man_phases:
                print("MAN files will be split into phase-based sweeps (upward/dwell/downward/turning)")
                print("Note: Phase splitting requires reading radar data to detect transitions")
            else:
                print("MAN files will be processed as multi-sweep netCDFs (one per aircraft track)")
            
            for f in man_files:
                try:
                    if split_man_phases:
                        # Read the file to analyze phase structure
                        print(f"Analyzing phases in: {os.path.basename(f)}")
                        radar_obj = read_mira35_mmclx(f, gzip_flag=gzip_flag, 
                                                     revised_northangle=revised_northangle)
                        
                        # Extract angle and time data
                        azimuth = radar_obj.azimuth['data']
                        elevation = radar_obj.elevation['data']
                        time_data = radar_obj.time['data']
                        
                        # Detect phase transitions
                        phases = detect_man_sweep_phases(azimuth, elevation, time_data)
                        print(f"  Detected {len(phases)} phase(s):")
                        for phase in phases:
                            n_rays = phase['end_idx'] - phase['start_idx'] + 1
                            print(f"    - {phase['phase']}: {n_rays} rays, "
                                  f"elev {phase['mean_el']:.1f}°, az {phase['mean_az']:.1f}°")
                        
                        # Apply phase-based sweep structure to radar object
                        radar_obj = apply_phase_sweeps_to_radar(radar_obj, phases)
                        print(f"  Applied {len(phases)} phase-based sweeps to radar object")
                        
                        # Write the modified radar object to file
                        # Create filename based on first ray time
                        from kepler_utils import cfradial_add_ncas_metadata, update_history_attribute
                        time_str = radar_obj.time['units'].split()[-1]
                        start_time = datetime.datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%SZ')
                        first_ray_time = start_time + datetime.timedelta(seconds=float(radar_obj.time['data'][0]))
                        
                        filename = f"ncas-mobile-ka-band-radar-1_{location}_{first_ray_time.strftime('%Y%m%d-%H%M%S')}_man_l1_v{data_version}.nc"
                        filepath = os.path.join(outdir, filename)
                        
                        # Verify radar object has required position fields
                        if not hasattr(radar_obj, 'latitude') or radar_obj.latitude is None:
                            print(f"  ERROR: Radar object missing latitude field")
                            continue
                        if not hasattr(radar_obj, 'longitude') or radar_obj.longitude is None:
                            print(f"  ERROR: Radar object missing longitude field")
                            continue
                        if not hasattr(radar_obj, 'altitude') or radar_obj.altitude is None:
                            print(f"  ERROR: Radar object missing altitude field")
                            continue
                        
                        # Confirm position data
                        lat = radar_obj.latitude['data'][0] if 'data' in radar_obj.latitude else radar_obj.latitude
                        lon = radar_obj.longitude['data'][0] if 'data' in radar_obj.longitude else radar_obj.longitude
                        alt = radar_obj.altitude['data'][0] if 'data' in radar_obj.altitude else radar_obj.altitude
                        print(f"  Position: lat={lat:.4f}, lon={lon:.4f}, alt={alt:.1f}m")
                        
                        # Write to netCDF
                        print(f"  Writing phase-split MAN file: {filename}")
                        print(f"  Radar dimensions: nrays={radar_obj.nrays}, ngates={radar_obj.ngates}, nsweeps={radar_obj.nsweeps}")
                        
                        # Debug: Check sweep dimensional arrays
                        print(f"  sweep_number shape: {radar_obj.sweep_number['data'].shape}")
                        print(f"  sweep_start_ray_index shape: {radar_obj.sweep_start_ray_index['data'].shape}")
                        print(f"  sweep_end_ray_index shape: {radar_obj.sweep_end_ray_index['data'].shape}")
                        print(f"  fixed_angle shape: {radar_obj.fixed_angle['data'].shape}")
                        print(f"  sweep_mode shape: {radar_obj.sweep_mode['data'].shape}")
                        if hasattr(radar_obj, 'target_scan_rate') and radar_obj.target_scan_rate is not None:
                            print(f"  target_scan_rate shape: {radar_obj.target_scan_rate['data'].shape}")
                        
                        # Debug: Check instrument parameter dimensions
                        if hasattr(radar_obj, 'instrument_parameters'):
                            for param_name, param_dict in radar_obj.instrument_parameters.items():
                                if 'data' in param_dict:
                                    print(f"  {param_name} shape: {param_dict['data'].shape}, dtype: {param_dict['data'].dtype}")
                        
                        try:
                            pyart.io.write_cfradial(filepath, radar_obj, format='NETCDF4')
                        except Exception as write_error:
                            print(f"  ERROR during write_cfradial: {write_error}")
                            import traceback
                            traceback.print_exc()
                            raise
                        
                        # Add elevation and azimuth scan rates for MAN scans
                        from kepler_utils import add_scan_rate_coordinates
                        add_scan_rate_coordinates(filepath, radar_obj)
                        
                        # Add NCAS metadata
                        print(f"  Adding NCAS metadata...")
                        print(f"    Project YAML: {yaml_project_file}")
                        print(f"    Instrument YAML: {yaml_instrument_file}")
                        print(f"    Tracking tag: {tracking_tag}")
                        print(f"    Data version: {data_version}")
                        
                        try:
                            cfradial_add_ncas_metadata(
                                filepath, 
                                yaml_project_file=yaml_project_file,
                                yaml_instrument_file=yaml_instrument_file,
                                tracking_tag=tracking_tag,
                                data_version=data_version
                            )
                            
                            # Add history entry with phase information
                            phase_summary = ', '.join([p['phase'] for p in phases])
                            history_msg = f"{campaign.upper()}: Phase-split MAN scan ({len(phases)} phases: {phase_summary})"
                            if revised_northangle:
                                history_msg += f", revised_northangle={revised_northangle}deg"
                            update_history_attribute(filepath, history_msg)
                            
                            print(f"  Successfully added metadata to {filename}")
                        except Exception as meta_error:
                            print(f"  ERROR: Failed to add NCAS metadata: {meta_error}")
                            import traceback
                            traceback.print_exc()
                            # Re-raise to make the error more visible
                            raise
                        
                    else:
                        # Standard processing: entire file as multi-sweep
                        RadarDS_MAN = multi_mmclx2cfrad(
                            [f], outdir, scan_name='MAN', gzip_flag=gzip_flag,
                            azimuth_offset=azimuth_offset,
                            tracking_tag=tracking_tag,
                            campaign=campaign,
                            revised_northangle=revised_northangle,
                            data_version=data_version,
                            single_sweep=False,  # Multi-sweep: captures all phases in one netCDF
                            yaml_project_file=yaml_project_file,
                            yaml_instrument_file=yaml_instrument_file
                        )
                    
                    print(f"Processed MAN file: {os.path.basename(f)} (aircraft tracking scan)")
                except Exception as e:
                    print(f"Error processing MAN file {f}: {e}")
                    import traceback
                    traceback.print_exc()
    except Exception as e:
        print(f"MAN processing problem: {e}")
        import traceback
        traceback.print_exc()
    
    print(f"Completed {campaign.upper()} processing for {datestr}")

def process_kepler_cobalt_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    azimuth_offset: float = 0.0,
    revised_northangle: float = 55.9,
    gzip_flag: bool = True,
    data_version: str = "1.0.3",
    single_sweep: bool = False,
    tracking_tag: str = 'AMOF_20231120125118',
    campaign: str = 'cobalt',  # Add campaign parameter
    no_vpt: bool = False  # Add no_vpt argument
) -> None:
    """
    Process COBALT campaign data for a single day - Step 1.
    """
    print(f"Processing {campaign.upper()} day: {datestr}")
    print(f"Single sweep mode: {single_sweep}")
    print(f"Data version: {data_version}")
    print(f"Tracking tag: {tracking_tag}")
    print(f"No VPT processing: {no_vpt}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Define the start and end times for the loop
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d')
    end_date = start_date + datetime.timedelta(days=1)
    
    print(f"start_date = {start_date} end_date = {end_date}")
    
    # Skip VPT processing if no_vpt is True
    if no_vpt:
        print("Skipping VPT file processing (--no-vpt option enabled)")
    else:
        print("Processing VPT files...")
        try:
            vpt_files = find_mmclxfiles(
                start_date.strftime('%Y-%m-%d %H:%M:%S'),
                end_date.strftime('%Y-%m-%d %H:%M:%S'),
                'vert', 
                inpath,
                gzip_flag=gzip_flag
            )
            print(vpt_files)
            print(f"Found {len(vpt_files)} VPT files")
            
            if len(vpt_files) > 0:
                print("VPT files will always create multi-sweep output (ignoring single_sweep flag)")
                RadarDS_VPT = multi_mmclx2cfrad(
                    vpt_files, outdir, scan_name='VPT', gzip_flag=gzip_flag,
                    azimuth_offset=azimuth_offset, 
                    tracking_tag=tracking_tag,
                    campaign=campaign,  # Use campaign parameter instead of hardcoded 'cobalt'
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,  # ALWAYS FALSE FOR VPT
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(vpt_files)} VPT files into multi-sweep dataset")
        except Exception as e:
            print(f"VPT problem: {e}")
            pass
    
    # Process RHI files (use single_sweep setting)
    print("Processing RHI files...")
    rhi_files = find_mmclx_rhi_files(
        start_date.strftime('%Y-%m-%d %H:%M:%S'), 
        end_date.strftime('%Y-%m-%d %H:%M:%S'), 
        0, 370, 
        inpath,
        gzip_flag=gzip_flag, 
        azimuth_offset=azimuth_offset,
        revised_northangle=revised_northangle
    )
    
    print(rhi_files)
    print(f"Found {len(rhi_files)} RHI files")
    
    if len(rhi_files) > 0:
        if single_sweep:
            # Process each RHI file separately
            for f in rhi_files:
                try:
                    RadarDS_RHI = multi_mmclx2cfrad(
                        [f], outdir, scan_name='RHI', gzip_flag=gzip_flag,
                        azimuth_offset=azimuth_offset, 
                        tracking_tag=tracking_tag,
                        campaign=campaign, 
                        revised_northangle=revised_northangle,
                        data_version=data_version,
                        single_sweep=True,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    print(f"Processed single RHI file: {f}")
                except Exception as e:
                    print(f"Error processing RHI file {f}: {e}")
        else:
            # Process all RHI files together (original behavior)
            try:
                RadarDS_RHI = multi_mmclx2cfrad(
                    rhi_files, outdir, scan_name='RHI', gzip_flag=gzip_flag,
                    azimuth_offset=azimuth_offset, 
                    tracking_tag=tracking_tag,
                    campaign=campaign, 
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(rhi_files)} RHI files into combined dataset")
            except Exception as e:
                print(f"Error processing RHI files: {e}")
    
    # Process VAD files for whole day - ALWAYS MULTI-SWEEP  
    print("Processing VAD files...")
    try:
        from kepler_utils import find_mmclx_vad_files
        vad_files = find_mmclx_vad_files(
            start_date.strftime('%Y-%m-%d %H:%M:%S'),
            end_date.strftime('%Y-%m-%d %H:%M:%S'),
            80, 90, 
            inpath,
            gzip_flag=gzip_flag
        )
        print(f"Found {len(vad_files)} VAD files")
        
        if len(vad_files) > 0:
            print("VAD files will always create multi-sweep output (ignoring single_sweep flag)")
            RadarDS_VAD = multi_mmclx2cfrad(
                vad_files, outdir, scan_name='VAD', gzip_flag=gzip_flag,
                azimuth_offset=azimuth_offset, 
                tracking_tag=tracking_tag,
                campaign=campaign, 
                revised_northangle=revised_northangle,
                data_version=data_version,
                single_sweep=False,  # ALWAYS FALSE FOR VAD
                yaml_project_file=yaml_project_file,
                yaml_instrument_file=yaml_instrument_file
            )
            print(f"Processed {len(vad_files)} VAD files into multi-sweep dataset")
    except Exception as e:
        print(f"VAD problem: {e}")
        pass
    
    print(f"Completed {campaign.upper()} processing for {datestr}")

def process_kepler_kasbex_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    azimuth_offset: float = 0.0,
    revised_northangle: float = 55.7,
    gzip_flag: bool = False,
    data_version: str = "1.0.0",
    single_sweep: bool = False,
    tracking_tag: str = 'AMOF_20250508133639',
    min_elevation: float = 0.0,
    max_elevation: float = 89.5,
    no_vpt: bool = False,  # Add no_vpt argument
    max_age: Optional[float] = None  # Add max_age argument
) -> None:
    """
    Process KASBEX campaign data for a single day - Step 1.
    Now handles VPT files with different range gate dimensions.
    """
    print(f"Processing KASBEX day: {datestr}")
    print(f"Using revised_northangle: {revised_northangle}° (from YAML configuration)")
    print(f"Single sweep mode: {single_sweep}")
    print(f"Data version: {data_version}")
    print(f"Tracking tag: {tracking_tag}")
    print(f"Input path: {inpath}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Define the start and end times for the loop
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d')
    end_date = start_date + datetime.timedelta(days=1)
    
    print(f"start_date = {start_date} end_date = {end_date}")
    
    # Check if inpath already includes the date (i.e., when using --latest)
    path_parts = inpath.rstrip('/').split('/')
    if datestr in path_parts or 'latest' in path_parts:
        # inpath already includes date/latest structure
        search_path = inpath
        print(f"Using direct path for file search: {search_path}")
        is_latest_mode = 'latest' in path_parts
    else:
        # Traditional structure - need to add date
        search_path = os.path.join(inpath, datestr)
        print(f"Using traditional date-based path: {search_path}")
        is_latest_mode = False
    
    # Process PPI files - prefer gzipped in latest mode
    print(f"Processing PPI files (elevation range: {min_elevation}° to {max_elevation}°)...")
    
    if is_latest_mode:
        print("Latest mode: Preferring gzipped PPI files")
        ppi_files = find_files_prefer_gzip(
        inpath, '*ppi*.mmclx*', force_gzip=gzip_flag, max_age_hours=max_age, use_max_age=True
        )
        #ppi_files = find_files_prefer_gzip(search_path, '*ppi*.mmclx*', force_gzip=True)
        # All files should be gzipped in this mode
        ppi_gzip_flag = True
    else:
        # Traditional mode - use glob and global gzip_flag
        import glob
        ppi_pattern = os.path.join(search_path, '*ppi*.mmclx*')
        ppi_files = glob.glob(ppi_pattern)
        ppi_files.sort()
        ppi_gzip_flag = gzip_flag
    
    print(f"Found {len(ppi_files)} PPI files")
    if ppi_files:
        print(f"Sample PPI files: {[os.path.basename(f) for f in ppi_files[:3]]}")
        if is_latest_mode:
            print(f"Using gzip mode for PPI files: {ppi_gzip_flag}")
    
    if len(ppi_files) > 0:
        if single_sweep:
            # Process each PPI file separately
            for f in ppi_files:
                try:
                    RadarDS_PPI = multi_mmclx2cfrad(
                        [f], outdir, scan_name='PPI', gzip_flag=ppi_gzip_flag,
                        azimuth_offset=azimuth_offset, 
                        tracking_tag=tracking_tag,
                        campaign='kasbex', 
                        revised_northangle=revised_northangle,
                        data_version=data_version,
                        single_sweep=True,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    print(f"Processed single PPI file: {os.path.basename(f)}")
                except Exception as e:
                    print(f"Error processing PPI file {f}: {e}")
        else:
            # Process all PPI files together
            try:
                RadarDS_PPI = multi_mmclx2cfrad(
                    ppi_files, outdir, scan_name='PPI', gzip_flag=ppi_gzip_flag,
                    azimuth_offset=azimuth_offset, 
                    tracking_tag=tracking_tag,
                    campaign='kasbex', 
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(ppi_files)} PPI files into combined dataset")
            except Exception as e:
                print(f"Error processing PPI files: {e}")
    
    # Process RHI files - prefer gzipped in latest mode
    print("Processing RHI files...")
    
    if is_latest_mode:
        print("Latest mode: Preferring gzipped RHI files")
        rhi_files = find_files_prefer_gzip(
            inpath, '*rhi*.mmclx*', force_gzip=gzip_flag, max_age_hours=max_age, use_max_age=True
        )
        print(f"Found {len(rhi_files)} RHI files")
        #rhi_files = find_files_prefer_gzip(search_path, '*rhi*.mmclx*', force_gzip=True)
        # All files should be gzipped in this mode
        rhi_gzip_flag = True
    else:
        # Traditional mode - use glob and global gzip_flag
        import glob
        rhi_pattern = os.path.join(search_path, '*rhi*.mmclx*')
        rhi_files = glob.glob(rhi_pattern)
        rhi_files.sort()
        rhi_gzip_flag = gzip_flag
    
    print(f"Found {len(rhi_files)} RHI files")
    if rhi_files:
        print(f"Sample RHI files: {[os.path.basename(f) for f in rhi_files[:3]]}")
        if is_latest_mode:
            print(f"Using gzip mode for RHI files: {rhi_gzip_flag}")
    
    if len(rhi_files) > 0:
        if single_sweep:
            # Process each RHI file separately
            for f in rhi_files:
                try:
                    RadarDS_RHI = multi_mmclx2cfrad(
                        [f], outdir, scan_name='RHI', gzip_flag=rhi_gzip_flag,
                        azimuth_offset=azimuth_offset, 
                        tracking_tag=tracking_tag,
                        campaign='kasbex', 
                        revised_northangle=revised_northangle,
                        data_version=data_version,
                        single_sweep=True,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    print(f"Processed single RHI file: {os.path.basename(f)}")
                except Exception as e:
                    print(f"Error processing RHI file {f}: {e}")
        else:
            # Process all RHI files together
            try:
                RadarDS_RHI = multi_mmclx2cfrad(
                    rhi_files, outdir, scan_name='RHI', gzip_flag=rhi_gzip_flag,
                    azimuth_offset=azimuth_offset, 
                    tracking_tag=tracking_tag,
                    campaign='kasbex', 
                    revised_northangle=revised_northangle,
                    data_version=data_version,
                    single_sweep=False,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                print(f"Processed {len(rhi_files)} RHI files into combined dataset")
            except Exception as e:
                print(f"Error processing RHI files: {e}")
    
    # Process VPT files with range gate grouping - allow both gzipped and uncompressed
    if no_vpt:
        print("Skipping VPT file processing (--no-vpt option enabled)")
    else:
        print("Processing VPT files...")
        try:
            if is_latest_mode:
                print("Latest mode: Allowing both gzipped and uncompressed VPT files")
                # Look for both compressed and uncompressed VPT files, but don't force gzip
                vpt_files = find_files_prefer_gzip(search_path, '*vert*.mmclx*', force_gzip=False)
            
                # Separate by compression status for processing
                vpt_files_gz = [f for f in vpt_files if f.endswith('.gz')]
                vpt_files_normal = [f for f in vpt_files if not f.endswith('.gz')]
            
                print(f"Found {len(vpt_files_gz)} gzipped VPT files")
                print(f"Found {len(vpt_files_normal)} normal VPT files")
                print(f"Total VPT files: {len(vpt_files)}")
            
                if vpt_files_gz:
                    print(f"Sample gzipped VPT files: {[os.path.basename(f) for f in vpt_files_gz[:3]]}")
                if vpt_files_normal:
                    print(f"Sample normal VPT files: {[os.path.basename(f) for f in vpt_files_normal[:3]]}")
            
                # Process gzipped VPT files with range gate grouping
                if len(vpt_files_gz) > 0:
                    print("Processing gzipped VPT files (multi-sweep output)")
                    vpt_files_gz.sort()
                
                    # Group by range gates
                    gz_file_groups = group_files_by_range_gates(vpt_files_gz, gzip_flag=True)
                
                    for n_gates, file_group in gz_file_groups.items():
                        if n_gates == 'error':
                            print(f"Skipping {len(file_group)} problematic gzipped VPT files")
                            continue
                        
                        print(f"Processing {len(file_group)} gzipped VPT files with {n_gates} range gates")
                        try:
                            RadarDS_VPT_gz = multi_mmclx2cfrad(
                                file_group, outdir, scan_name='VPT', gzip_flag=True,
                                azimuth_offset=azimuth_offset, 
                                tracking_tag=tracking_tag,
                                campaign='kasbex', 
                                revised_northangle=revised_northangle,
                                data_version=data_version,
                                single_sweep=False,  # Always FALSE for VPT
                                yaml_project_file=yaml_project_file,
                                yaml_instrument_file=yaml_instrument_file
                            )
                            print(f"Successfully processed {len(file_group)} gzipped VPT files ({n_gates} gates)")
                        except Exception as e:
                            print(f"Error processing gzipped VPT files with {n_gates} gates: {e}")
            
                # Process normal VPT files with range gate grouping
                if len(vpt_files_normal) > 0:
                    print("Processing uncompressed VPT files (multi-sweep output)")
                    vpt_files_normal.sort()
                
                    # Group by range gates
                    normal_file_groups = group_files_by_range_gates(vpt_files_normal, gzip_flag=False)
                
                    for n_gates, file_group in normal_file_groups.items():
                        if n_gates == 'error':
                            print(f"Skipping {len(file_group)} problematic normal VPT files")
                            continue
                        
                        print(f"Processing {len(file_group)} normal VPT files with {n_gates} range gates")
                        try:
                            RadarDS_VPT_normal = multi_mmclx2cfrad(
                                file_group, outdir, scan_name='VPT', gzip_flag=False,
                                azimuth_offset=azimuth_offset, 
                                tracking_tag=tracking_tag,
                                campaign='kasbex', 
                                revised_northangle=revised_northangle,
                                data_version=data_version,
                                single_sweep=False,  # Always FALSE for VPT
                                yaml_project_file=yaml_project_file,
                                yaml_instrument_file=yaml_instrument_file
                            )
                            print(f"Successfully processed {len(file_group)} normal VPT files ({n_gates} gates)")
                        except Exception as e:
                            print(f"Error processing normal VPT files with {n_gates} gates: {e}")
        
            else:
                # Traditional mode - use global gzip_flag for all VPT files with range gate grouping
                import glob
                vpt_pattern = os.path.join(search_path, '*vert*.mmclx*')
                vpt_files = glob.glob(vpt_pattern)
                vpt_files.sort()
            
                print(f"Found {len(vpt_files)} VPT files")
                if vpt_files:
                    print(f"Sample VPT files: {[os.path.basename(f) for f in vpt_files[:3]]}")
            
                if len(vpt_files) > 0:
                    print("Processing VPT files (multi-sweep output) with range gate grouping")
                
                    # Group by range gates
                    file_groups = group_files_by_range_gates(vpt_files, gzip_flag=gzip_flag)
                
                    for n_gates, file_group in file_groups.items():
                        if n_gates == 'error':
                            print(f"Skipping {len(file_group)} problematic VPT files")
                            continue
                        
                        print(f"Processing {len(file_group)} VPT files with {n_gates} range gates")
                        try:
                            RadarDS_VPT = multi_mmclx2cfrad(
                                file_group, outdir, scan_name='VPT', gzip_flag=gzip_flag,
                                azimuth_offset=azimuth_offset, 
                                tracking_tag=tracking_tag,
                                campaign='kasbex', 
                                revised_northangle=revised_northangle,
                                data_version=data_version,
                                single_sweep=False,  # Always FALSE for VPT
                                yaml_project_file=yaml_project_file,
                                yaml_instrument_file=yaml_instrument_file
                            )
                            print(f"Successfully processed {len(file_group)} VPT files ({n_gates} gates)")
                        except Exception as e:
                            print(f"Error processing VPT files with {n_gates} gates: {e}")
            
        except Exception as e:
            print(f"VPT problem: {e}")
            import traceback
            traceback.print_exc()
            pass
    
    print(f"Completed KASBEX processing for {datestr} with north_angle={revised_northangle}°")

def find_files_prefer_gzip(search_path, pattern, force_gzip=False, max_age_hours=None, use_max_age=False):
    """
    Find files matching pattern, preferring gzipped versions over uncompressed ones,
    and optionally filter by maximum age using timestamps from filenames.

    Args:
        search_path (str): Directory to search in
        pattern (str): Glob pattern for files (e.g., '*ppi*.mmclx*')
        force_gzip (bool): If True, only return gzipped files
        max_age_hours (float): Maximum file age in hours (None for no limit)
        use_max_age (bool): Whether to apply max-age filtering (only if --latest is used)

    Returns:
        list: List of files with gzipped versions preferred and filtered by age
    """
    import glob
    import os
    import datetime
    import re

    # Find all matching files
    all_files = glob.glob(os.path.join(search_path, pattern))
    all_files.sort()

    # Apply max-age filtering only if use_max_age is True
    if use_max_age and max_age_hours is not None:
        current_time = datetime.datetime.utcnow()
        max_age_seconds = max_age_hours * 3600
        filtered_files = []
        time_pattern = r'(\d{8})-(\d{6})'  # Match YYYYMMDD-HHMMSS in filenames

        for file_path in all_files:
            try:
                filename = os.path.basename(file_path)
                match = re.search(time_pattern, filename)
                if match:
                    date_str = match.group(1)  # YYYYMMDD
                    time_str = match.group(2)  # HHMMSS
                    file_time = datetime.datetime.strptime(f"{date_str}-{time_str}", "%Y%m%d-%H%M%S")
                    file_age_seconds = (current_time - file_time).total_seconds()

                    if file_age_seconds <= max_age_seconds:
                        filtered_files.append(file_path)
                    else:
                        print(f"Skipping old file: {filename} (age: {file_age_seconds / 3600:.1f} hours)")
                else:
                    print(f"Skipping file with no timestamp: {filename}")
            except Exception as e:
                print(f"Error checking file age for {file_path}: {e}")
                # Include file if age check fails
                filtered_files.append(file_path)
        all_files = filtered_files

    if force_gzip:
        # Only return gzipped files
        gzip_files = [f for f in all_files if f.endswith('.gz')]
        print(f"Force gzip mode: found {len(gzip_files)} gzipped files out of {len(all_files)} total")
        return gzip_files

    # Create dictionary to group files by base name
    file_groups = {}

    for file_path in all_files:
        if file_path.endswith('.gz'):
            # Remove .gz to get base name
            base_name = file_path[:-3]
            if base_name not in file_groups:
                file_groups[base_name] = {'gz': None, 'normal': None}
            file_groups[base_name]['gz'] = file_path
        else:
            # This is the base name
            base_name = file_path
            if base_name not in file_groups:
                file_groups[base_name] = {'gz': None, 'normal': None}
            file_groups[base_name]['normal'] = file_path

    # Build final file list, preferring gzipped versions
    final_files = []
    for base_name, versions in file_groups.items():
        if versions['gz'] is not None:
            # Prefer gzipped version
            final_files.append(versions['gz'])
            if versions['normal'] is not None:
                print(f"Found both versions, preferring gzipped: {os.path.basename(versions['gz'])}")
        elif versions['normal'] is not None:
            # Use uncompressed version
            final_files.append(versions['normal'])

    final_files.sort()
    return final_files

def load_yaml_config(yaml_project_file, yaml_instrument_file):
    """
    Load configuration from YAML files, including north_angle.
    
    Args:
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        
    Returns:
        Dictionary of configuration parameters. north_angle is None if not found in YAML.
    """
    config = {
        'north_angle': None,  # Will be None if not found in YAML
        'data_version': '1.0.0',
        'location': 'cao'
    }
    
    # Load project YAML
    try:
        with open(yaml_project_file, 'r') as f:
            project_data = yaml.safe_load(f)
        
        print(f"Loading configuration from: {yaml_project_file}")
        
        # Navigate through your YAML structure: ncas_instruments -> ncas-mobile-ka-band-radar-1
        if 'ncas_instruments' in project_data:
            for instrument in project_data['ncas_instruments']:
                if 'ncas-mobile-ka-band-radar-1' in instrument:
                    radar_config = instrument['ncas-mobile-ka-band-radar-1']
                    
                    # Extract north_angle
                    if 'north_angle' in radar_config:
                        config['north_angle'] = float(radar_config['north_angle'])
                        print(f"Loaded north_angle from YAML: {config['north_angle']}°")
                    
                    # Extract other parameters if available
                    if 'processing_software' in radar_config and 'version' in radar_config['processing_software']:
                        config['data_version'] = radar_config['processing_software']['version']
                    
                    if 'platform' in radar_config and 'location' in radar_config['platform']:
                        config['location'] = radar_config['platform']['location'].lower()
                    
                    break
        
    except Exception as e:
        print(f"Warning: Could not load project YAML config: {e}")
    
    return config

def get_campaign_info(campaign: str) -> Dict[str, Any]:
    """
    Get campaign-specific configuration parameters.
    
    Args:
        campaign: Campaign name
        
    Returns:
        Dictionary of campaign-specific parameters
    """
    campaign_configs = {
        'cobalt': {
            'revised_northangle': 55.62,
            'tracking_tag': 'AMOF_20231120125118'
        },
        'woest': {
            'revised_northangle': 302.15,  # Different north angle for WOEST
            'tracking_tag': 'AMOF_20220922221548'
        },
        'ccrest-m': {
            'revised_northangle': 55.7,
            'tracking_tag': 'AMOF_20230201132601'
        },
        'coalesc3': {
            'revised_northangle': 287.0,
            'tracking_tag': 'AMF_07092016101810'
        },
        'kasbex': {
            'revised_northangle': 55.62,
            'tracking_tag': 'AMOF_20250508133639'  # KASBEX-specific tracking tag
        },
        'picasso': {
            'revised_northangle': 97.1,  # PICASSO north angle from YAML
            'tracking_tag': 'CFARR_0002'
        }
    }
    
    return campaign_configs.get(campaign, {
        'revised_northangle': 55.62,  # Default
        'tracking_tag': f'AMOF_{campaign.upper()}'
    })

def group_files_by_range_gates(files, gzip_flag=False):
    """
    Group mmclx files by their range gate dimensions to avoid concatenation errors.
    
    Args:
        files (list): List of mmclx file paths
        gzip_flag (bool): Whether files are gzipped
        
    Returns:
        dict: Dictionary with range_gates as keys and lists of files as values
    """
    from kepler_utils import read_mira35_mmclx
    import os
    
    file_groups = {}
    
    for file_path in files:
        try:
            # Read just the first sweep to get dimensions
            radar_ds = read_mira35_mmclx(file_path, gzip_flag=gzip_flag, revised_northangle=55.62)
            n_gates = radar_ds.ngates
            
            if n_gates not in file_groups:
                file_groups[n_gates] = []
            file_groups[n_gates].append(file_path)
            
            print(f"File {os.path.basename(file_path)}: {n_gates} gates")
            
        except Exception as e:
            print(f"Warning: Could not read {file_path} for grouping: {e}")
            # Put problematic files in a separate group
            if 'error' not in file_groups:
                file_groups['error'] = []
            file_groups['error'].append(file_path)
    
    return file_groups

# Example usage functions
def process_campaign_day(
    campaign: str,
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    gzip_flag: bool = True,
    data_version: str = "1.0.0",
    single_sweep: bool = False,
    revised_northangle: float = None,
    no_vpt: bool = False,  # Add no_vpt argument
    max_age: Optional[float] = None,
    **kwargs
) -> None:
    """
    Process a single day of campaign data.
    
    Args:
        campaign: Campaign name ('woest', 'ccrest-m', 'coalesc3', 'picasso', 'cobalt', 'kasbex')
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        outpath: Output directory path
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        gzip_flag: Whether input files are gzip compressed
        data_version: Data version string
        single_sweep: If True, create separate files for each sweep
        revised_northangle: North angle override (if None, load from YAML, or use campaign default)
        **kwargs: Additional campaign-specific arguments
    """
    
    # Get campaign-specific configuration first (provides defaults)
    campaign_lower = campaign.lower()
    campaign_info = get_campaign_info(campaign_lower)
    
    # Load configuration from YAML
    yaml_config = load_yaml_config(yaml_project_file, yaml_instrument_file)
    
    # Determine north angle with proper priority:
    # 1. Command line override (if provided)
    # 2. YAML config (if north_angle exists in YAML)
    # 3. Campaign-specific default from campaign_info
    if revised_northangle is not None:
        # Command line override
        print(f"Using north_angle from command line: {revised_northangle}°")
    elif yaml_config.get('north_angle') is not None:
        # YAML has a specific value
        revised_northangle = yaml_config['north_angle']
        print(f"Using north_angle from YAML: {revised_northangle}°")
    else:
        # Use campaign-specific default
        revised_northangle = campaign_info.get('revised_northangle', 55.62)
        print(f"Using north_angle from campaign default: {revised_northangle}°")
    
    # Define available processors
    processors = {
        'woest': process_kepler_woest_day_step1,
        'ccrest-m': process_kepler_general_day_step1,
        'coalesc3': process_kepler_general_day_step1,
        'picasso': process_kepler_picasso_day_step1,
        'cobalt': process_kepler_cobalt_day_step1,
        'kasbex': process_kepler_kasbex_day_step1
    }
    
    if campaign_lower not in processors:
        raise ValueError(f"Unknown campaign: {campaign}. Available: {list(processors.keys())}")
    
    # Set up processing arguments
    processing_args = {
        'datestr': datestr,
        'inpath': inpath,
        'outpath': outpath,
        'yaml_project_file': yaml_project_file,
        'yaml_instrument_file': yaml_instrument_file,
        'gzip_flag': gzip_flag,
        'data_version': data_version,
        'single_sweep': single_sweep,
        'revised_northangle': revised_northangle,  # Use the north angle from YAML or override
        'tracking_tag': campaign_info.get('tracking_tag', f'AMOF_{campaign.upper()}'),
        'campaign': campaign_lower,  # Pass campaign name for general processor
        'no_vpt': no_vpt,  # Pass no_vpt argument
        **kwargs
    }
    
    print(f"Processing {campaign.upper()} campaign data for {datestr}")
    print(f"Configuration:")
    print(f"  North angle: {revised_northangle}°")
    print(f"  Data version: {data_version}")
    print(f"  Single sweep: {single_sweep}")
    print(f"  Tracking tag: {processing_args['tracking_tag']}")
    
    # Call the appropriate processor
    processor = processors[campaign_lower]
    processor(**processing_args)
    
    print(f"Completed {campaign.upper()} processing for {datestr}")