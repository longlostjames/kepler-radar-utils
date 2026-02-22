#!/bin/bash
# cleanup_kepler_kasbex_files.sh
# Script to manage KASBEX processed files - move previous versions to "_old" and clean up numbered files

# Default values
DATE=""
BASE_PATH="/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/kasbex/L1c"
DRY_RUN=false
VERBOSE=false

# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Manage KASBEX processed files - move previous versions to '_old' and clean up numbered files"
    echo ""
    echo "OPTIONS:"
    echo "  -d, --date DATE        Date in YYYYMMDD format (default: today)"
    echo "  -p, --path PATH        Base path to KASBEX L1c data (default: $BASE_PATH)"
    echo "  --dry-run              Show what would be done without making changes"
    echo "  -v, --verbose          Show detailed output"
    echo "  -h, --help             Show this help message"
    echo ""
    echo "EXAMPLES:"
    echo "  $0 -d 20250801                    # Process files for specific date"
    echo "  $0 --dry-run                      # Show what would be done for today"
    echo "  $0 -d 20250801 --verbose          # Process with detailed output"
    echo ""
}

# Function to log messages
log_message() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case $level in
        "INFO")
            echo "[$timestamp] INFO: $message"
            ;;
        "WARN")
            echo "[$timestamp] WARN: $message" >&2
            ;;
        "ERROR")
            echo "[$timestamp] ERROR: $message" >&2
            ;;
        "DEBUG")
            if [ "$VERBOSE" = true ]; then
                echo "[$timestamp] DEBUG: $message"
            fi
            ;;
    esac
}

# Function to move file with backup
safe_move() {
    local src=$1
    local dst=$2
    
    if [ "$DRY_RUN" = true ]; then
        echo "DRY RUN: Would move '$src' -> '$dst'"
        return 0
    fi
    
    if [ -f "$src" ]; then
        if [ -f "$dst" ]; then
            log_message "WARN" "Destination file already exists: $dst"
            # Create backup with timestamp
            local backup="${dst}.backup.$(date +%Y%m%d_%H%M%S)"
            log_message "INFO" "Creating backup: $backup"
            mv "$dst" "$backup"
        fi
        
        log_message "INFO" "Moving: $(basename $src) -> $(basename $dst)"
        mv "$src" "$dst"
        
        if [ $? -eq 0 ]; then
            log_message "DEBUG" "Successfully moved $src to $dst"
            return 0
        else
            log_message "ERROR" "Failed to move $src to $dst"
            return 1
        fi
    else
        log_message "WARN" "Source file does not exist: $src"
        return 1
    fi
}

# Function to process files in a directory
process_date_directory() {
    local date_dir=$1
    local total_processed=0
    local total_errors=0
    
    if [ ! -d "$date_dir" ]; then
        log_message "ERROR" "Directory does not exist: $date_dir"
        return 1
    fi
    
    log_message "INFO" "Processing directory: $date_dir"
    
    # Change to the date directory
    cd "$date_dir" || {
        log_message "ERROR" "Cannot change to directory: $date_dir"
        return 1
    }
    
    # Find all numbered files (ending with _01, _02, etc.)
    local numbered_files=$(ls *_[0-9][0-9].nc 2>/dev/null)
    
    if [ -z "$numbered_files" ]; then
        log_message "INFO" "No numbered files found in $date_dir"
        return 0
    fi
    
    log_message "INFO" "Found numbered files:"
    for file in $numbered_files; do
        log_message "DEBUG" "  - $file"
    done
    
    # Group files by base name (everything before the _XX.nc)
    declare -A file_groups
    
    for file in $numbered_files; do
        # Extract base name by removing the _XX.nc suffix
        base_name=$(echo "$file" | sed 's/_[0-9][0-9]\.nc$/.nc/')
        
        if [ -z "${file_groups[$base_name]}" ]; then
            file_groups[$base_name]="$file"
        else
            file_groups[$base_name]="${file_groups[$base_name]} $file"
        fi
    done
    
    # Process each group
    for base_name in "${!file_groups[@]}"; do
        log_message "INFO" "Processing file group: $base_name"
        
        # Get all files for this base name, sorted
        files=($(echo ${file_groups[$base_name]} | tr ' ' '\n' | sort))
        
        log_message "DEBUG" "Files in group: ${files[*]}"
        
        # If base file already exists, move it to _old version
        if [ -f "$base_name" ]; then
            old_name=$(echo "$base_name" | sed 's/\.nc$/_old.nc/')
            log_message "INFO" "Moving existing file to old version"
            
            if safe_move "$base_name" "$old_name"; then
                ((total_processed++))
            else
                ((total_errors++))
                continue
            fi
        fi
        
        # Get the highest numbered file (should be the latest)
        latest_file="${files[-1]}"
        log_message "INFO" "Latest file identified: $latest_file"
        
        # Move the latest numbered file to the base name
        if safe_move "$latest_file" "$base_name"; then
            ((total_processed++))
        else
            ((total_errors++))
            continue
        fi
        
        # Remove other numbered files (they're older versions)
        for file in "${files[@]}"; do
            if [ "$file" != "$latest_file" ]; then
                if [ "$DRY_RUN" = true ]; then
                    echo "DRY RUN: Would remove old numbered file: $file"
                else
                    log_message "INFO" "Removing old numbered file: $file"
                    rm -f "$file"
                    if [ $? -eq 0 ]; then
                        ((total_processed++))
                    else
                        log_message "ERROR" "Failed to remove $file"
                        ((total_errors++))
                    fi
                fi
            fi
        done
    done
    
    log_message "INFO" "Completed processing $date_dir"
    log_message "INFO" "Total operations: $total_processed successful, $total_errors errors"
    
    return $total_errors
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--date)
            DATE="$2"
            shift 2
            ;;
        -p|--path)
            BASE_PATH="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

# Set default date to today if not provided
if [ -z "$DATE" ]; then
    DATE=$(date +%Y%m%d)
fi

# Validate date format
if ! echo "$DATE" | grep -qE '^[0-9]{8}$'; then
    log_message "ERROR" "Invalid date format. Use YYYYMMDD. Got: $DATE"
    exit 1
fi

# Validate date
if ! date -d "$DATE" >/dev/null 2>&1; then
    log_message "ERROR" "Invalid date: $DATE"
    exit 1
fi

# Check if base path exists
if [ ! -d "$BASE_PATH" ]; then
    log_message "ERROR" "Base path does not exist: $BASE_PATH"
    exit 1
fi

# Construct full date directory path
DATE_DIR="$BASE_PATH/$DATE"

log_message "INFO" "KASBEX File Cleanup Script Started"
log_message "INFO" "Date: $DATE"
log_message "INFO" "Base path: $BASE_PATH"
log_message "INFO" "Date directory: $DATE_DIR"
log_message "INFO" "Dry run: $DRY_RUN"
log_message "INFO" "Verbose: $VERBOSE"

# Process the date directory
if process_date_directory "$DATE_DIR"; then
    log_message "INFO" "File cleanup completed successfully"
    exit 0
else
    log_message "ERROR" "File cleanup completed with errors"
    exit 1
fi