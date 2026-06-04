#!/usr/bin/env python3
"""
generate_quicklooks_index.py

Generates a browsable index.html for a Kepler radar quicklooks directory.
Run from anywhere; pass -o to specify the quicklooks root.

Usage:
    python generate_quicklooks_index.py -o /data/processing/kepler/reading-general/quicklooks
"""

import argparse
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


# ---------------------------------------------------------------------------
# Filename parsing
# ---------------------------------------------------------------------------

FILENAME_RE = re.compile(
    r"(?P<instrument>[^_]+(?:_[^_]+)*)_(?P<date>\d{8})(?:-\d{6})?_(?P<scantype>[a-z][a-z0-9-]*)_"
    r"(?:[^_]+_)*(?P<level>l\d+)_(?P<version>v[\d.]+)\.png$",
    re.IGNORECASE,
)

SCAN_ORDER = ["vpt", "ppi", "ppi-map", "rhi"]


def parse_filename(fname):
    """Return (date_str, scantype) or None if unrecognised."""
    m = FILENAME_RE.search(fname)
    if m:
        return m.group("date"), m.group("scantype").lower()
    return None


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

def discover_images(root: Path):
    """Return dict: {date_str: {scantype: [rel_path, ...]}}."""
    data = defaultdict(lambda: defaultdict(list))
    for png in sorted(root.rglob("*.png")):
        rel = png.relative_to(root)
        parsed = parse_filename(png.name)
        if parsed:
            date_str, scantype = parsed
            data[date_str][scantype].append(str(rel))
    return data


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------

HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<style>
  :root {{
    --bg: #1a1a2e;
    --surface: #16213e;
    --card: #0f3460;
    --accent: #e94560;
    --text: #eaeaea;
    --muted: #8892a4;
    --tab-active: #e94560;
    --tab-inactive: #0f3460;
    --border: #253a5e;
    --sidebar-w: 200px;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: 'Segoe UI', system-ui, sans-serif;
    background: var(--bg);
    color: var(--text);
    min-height: 100vh;
  }}
  header {{
    background: var(--surface);
    border-bottom: 2px solid var(--accent);
    padding: 1rem 2rem;
    display: flex;
    align-items: center;
    gap: 1rem;
    position: sticky;
    top: 0;
    z-index: 100;
  }}
  header h1 {{ font-size: 1.4rem; font-weight: 600; }}
  header .sub {{ color: var(--muted); font-size: 0.9rem; }}
  .generated {{ margin-left: auto; color: var(--muted); font-size: 0.8rem; }}

  .layout {{
    display: flex;
    align-items: flex-start;
  }}

  /* Sidebar */
  .sidebar {{
    width: var(--sidebar-w);
    flex-shrink: 0;
    position: sticky;
    top: 60px;
    max-height: calc(100vh - 60px);
    overflow-y: auto;
    padding: 1rem 0.75rem;
    background: var(--surface);
    border-right: 1px solid var(--border);
  }}
  .sidebar-title {{
    font-size: 0.7rem;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    color: var(--muted);
    margin-bottom: 0.6rem;
  }}
  .year-group {{ margin-bottom: 0.4rem; }}
  .year-toggle {{
    width: 100%;
    background: none;
    border: none;
    color: var(--text);
    font-size: 0.88rem;
    font-weight: 600;
    text-align: left;
    cursor: pointer;
    padding: 0.25rem 0.3rem;
    border-radius: 3px;
    display: flex;
    justify-content: space-between;
    align-items: center;
  }}
  .year-toggle:hover {{ background: var(--card); }}
  .year-toggle .arrow {{ font-size: 0.65rem; transition: transform 0.2s; }}
  .year-toggle.open .arrow {{ transform: rotate(90deg); }}
  .month-list {{ display: none; padding-left: 0.6rem; }}
  .month-list.open {{ display: block; }}
  .month-link {{
    display: block;
    padding: 0.2rem 0.3rem;
    font-size: 0.82rem;
    color: var(--muted);
    text-decoration: none;
    border-radius: 3px;
  }}
  .month-link:hover {{ background: var(--card); color: var(--text); }}
  .month-link.active {{ color: var(--accent); font-weight: 600; }}

  main {{ padding: 1.5rem 2rem; flex: 1; min-width: 0; }}

  /* Scan-type tab bar */
  .tab-bar {{
    display: flex;
    gap: 0.5rem;
    margin-bottom: 1.5rem;
    flex-wrap: wrap;
  }}
  .tab-btn {{
    padding: 0.4rem 1.2rem;
    border: 1px solid var(--border);
    border-radius: 4px;
    background: var(--tab-inactive);
    color: var(--text);
    cursor: pointer;
    font-size: 0.9rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    transition: background 0.15s;
  }}
  .tab-btn.active {{
    background: var(--tab-active);
    border-color: var(--tab-active);
    font-weight: 600;
  }}
  .tab-btn:hover:not(.active) {{ background: var(--card); }}

  /* Year / month hierarchy */
  .year-section {{ margin-bottom: 2.5rem; }}
  .year-heading {{
    font-size: 1.15rem;
    font-weight: 700;
    color: var(--accent);
    border-bottom: 2px solid var(--border);
    padding-bottom: 0.4rem;
    margin-bottom: 1.2rem;
    letter-spacing: 0.04em;
    cursor: pointer;
    user-select: none;
  }}
  .year-heading .toggle-arrow {{ font-size: 0.75rem; margin-right: 0.5rem; transition: transform 0.2s; }}
  .year-heading.collapsed .toggle-arrow {{ transform: rotate(-90deg); }}
  .year-content.collapsed {{ display: none; }}
  .month-section {{ margin-bottom: 1.5rem; }}
  .month-heading {{
    font-size: 1rem;
    font-weight: 600;
    color: var(--text);
    border-bottom: 1px solid var(--border);
    padding-bottom: 0.3rem;
    margin-bottom: 0.8rem;
    letter-spacing: 0.03em;
    cursor: pointer;
    user-select: none;
  }}
  .month-heading .toggle-arrow {{ font-size: 0.75rem; margin-right: 0.5rem; transition: transform 0.2s; }}
  .month-heading.collapsed .toggle-arrow {{ transform: rotate(-90deg); }}
  .month-content.collapsed {{ display: none; }}
  /* Date sections */
  .date-section {{ margin-bottom: 1.5rem; }}
  .date-heading {{
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--muted);
    padding-bottom: 0.3rem;
    margin-bottom: 0.6rem;
    letter-spacing: 0.04em;
  }}
  .date-heading.collapsible {{
    cursor: pointer;
    user-select: none;
  }}
  .date-heading .toggle-arrow {{ font-size: 0.65rem; margin-right: 0.4rem; transition: transform 0.2s; }}
  .date-heading.collapsed .toggle-arrow {{ transform: rotate(-90deg); }}
  .date-content.collapsed {{ display: none; }}
  .image-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(340px, 1fr));
    gap: 1rem;
  }}
  .image-card {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 6px;
    overflow: hidden;
    transition: transform 0.15s, box-shadow 0.15s;
  }}
  .image-card:hover {{
    transform: translateY(-2px);
    box-shadow: 0 4px 18px rgba(0,0,0,0.4);
  }}
  .image-card a {{ display: block; }}
  .image-card img {{
    width: 100%;
    height: auto;
    display: block;
  }}
  .image-caption {{
    padding: 0.5rem 0.75rem;
    font-size: 0.75rem;
    color: var(--muted);
    word-break: break-all;
  }}

  /* Scan-type panels */
  .scan-panel {{ display: none; }}
  .scan-panel.active {{ display: block; }}

  /* Empty state */
  .empty {{ color: var(--muted); padding: 2rem 0; }}
</style>
</head>
<body>
<header>
  <div>
    <h1>{title}</h1>
    <div class="sub">{subtitle}</div>
  </div>
  <div class="generated">Generated {generated}</div>
</header>
<div class="layout">
<nav class="sidebar" id="sidebar">
  <div class="sidebar-title">Browse</div>
{sidebar_nav}
</nav>
<main>

<div class="tab-bar" id="tab-bar">
{tab_buttons}
</div>

{scan_panels}

</main>
</div>
<script>
  // Tab switching
  const tabBar = document.getElementById('tab-bar');
  const panels = document.querySelectorAll('.scan-panel');
  function updateSidebarHrefs(scan) {{
    document.querySelectorAll('.month-link').forEach(l => {{
      l.setAttribute('href', l.getAttribute('href').replace(/#panel-[^-]+-/, '#panel-' + scan + '-'));
    }});
  }}
  tabBar.addEventListener('click', e => {{
    const btn = e.target.closest('.tab-btn');
    if (!btn) return;
    const target = btn.dataset.target;
    const scan = target.replace('panel-', '');
    tabBar.querySelectorAll('.tab-btn').forEach(b => b.classList.toggle('active', b === btn));
    panels.forEach(p => p.classList.toggle('active', p.id === target));
    updateSidebarHrefs(scan);
    updateSidebarLinks();
  }});

  // Year collapse in main content
  document.querySelectorAll('.year-heading').forEach(h => {{
    h.addEventListener('click', () => {{
      h.classList.toggle('collapsed');
      h.nextElementSibling.classList.toggle('collapsed');
    }});
  }});

  // Month collapse in main content — collapse all except the latest in each panel
  document.querySelectorAll('.scan-panel').forEach(panel => {{
    const months = panel.querySelectorAll('.month-heading');
    months.forEach((h, i) => {{
      if (i > 0) {{
        h.classList.add('collapsed');
        h.nextElementSibling.classList.add('collapsed');
      }}
    }});
  }});
  document.querySelectorAll('.month-heading').forEach(h => {{
    h.addEventListener('click', () => {{
      h.classList.toggle('collapsed');
      h.nextElementSibling.classList.toggle('collapsed');
    }});
  }});

  // Day collapse in main content — collapse all except the latest in each panel
  document.querySelectorAll('.scan-panel').forEach(panel => {{
    const days = panel.querySelectorAll('.date-heading.collapsible');
    days.forEach((h, i) => {{
      if (i > 0) {{
        h.classList.add('collapsed');
        h.nextElementSibling.classList.add('collapsed');
      }}
    }});
  }});
  document.querySelectorAll('.date-heading.collapsible').forEach(h => {{
    h.addEventListener('click', () => {{
      h.classList.toggle('collapsed');
      h.nextElementSibling.classList.toggle('collapsed');
    }});
  }});

  // Sidebar year toggles
  document.querySelectorAll('.year-toggle').forEach(btn => {{
    btn.addEventListener('click', () => {{
      btn.classList.toggle('open');
      btn.nextElementSibling.classList.toggle('open');
    }});
  }});

  // Highlight active sidebar month link based on scroll position
  function updateSidebarLinks() {{
    const activePanel = document.querySelector('.scan-panel.active');
    if (!activePanel) return;
    const sections = activePanel.querySelectorAll('.month-section[data-month-id]');
    const links = document.querySelectorAll('.month-link');
    links.forEach(l => l.classList.remove('active'));
    let current = null;
    sections.forEach(s => {{
      if (s.getBoundingClientRect().top <= 120) current = s.dataset.monthId;
    }});
    if (current) {{
      const link = document.querySelector(`.month-link[href="#${{current}}"]`);
      if (link) link.classList.add('active');
    }}
  }}
  window.addEventListener('scroll', updateSidebarLinks, {{ passive: true }});
  updateSidebarLinks();

  // Open the first year in the sidebar by default
  const firstToggle = document.querySelector('.year-toggle');
  if (firstToggle) {{ firstToggle.classList.add('open'); firstToggle.nextElementSibling.classList.add('open'); }}
</script>
</body>
</html>
"""


MONTH_NAMES = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December",
]


def format_date(date_str):
    try:
        return datetime.strptime(date_str, "%Y%m%d").strftime("%-d %B %Y")
    except ValueError:
        return date_str


def group_dates_by_year_month(all_dates):
    """Return OrderedDict: {year_str: {month_str: [date_str, ...]}} newest-first."""
    from collections import OrderedDict
    ym = defaultdict(lambda: defaultdict(list))
    for d in all_dates:
        ym[d[:4]][d[4:6]].append(d)
    result = OrderedDict()
    for year in sorted(ym, reverse=True):
        result[year] = OrderedDict()
        for month in sorted(ym[year], reverse=True):
            result[year][month] = sorted(ym[year][month], reverse=True)
    return result


def build_scan_panel(scantype, by_scan, year_month_dates):
    """Build a <div class="scan-panel"> grouped by year then month."""
    year_sections = []
    for year, months in year_month_dates.items():
        month_sections = []
        for month, dates in months.items():
            month_name = MONTH_NAMES[int(month) - 1]
            month_id = f"panel-{scantype}-{year}-{month}"
            date_sections = []
            for date_str in dates:
                images = by_scan.get(scantype, {}).get(date_str, [])
                if not images:
                    continue
                cards = "".join(
                    f'      <div class="image-card">'
                    f'<a href="{rel_path}" target="_blank">'
                    f'<img src="{rel_path}" alt="{os.path.basename(rel_path)}" loading="lazy"></a>'
                    f'<div class="image-caption">{os.path.basename(rel_path)}</div>'
                    f"</div>\n"
                    for rel_path in images
                )
                if scantype in ('ppi', 'ppi-map', 'rhi'):
                    date_sections.append(
                        f'    <div class="date-section">\n'
                        f'      <div class="date-heading collapsible"><span class="toggle-arrow">▼</span>{format_date(date_str)}</div>\n'
                        f'      <div class="date-content">\n'
                        f'      <div class="image-grid">\n{cards}      </div>\n'
                        f'      </div>\n'
                        f"    </div>"
                    )
                else:
                    date_sections.append(
                        f'    <div class="date-section">\n'
                        f'      <div class="date-heading">{format_date(date_str)}</div>\n'
                        f'      <div class="image-grid">\n{cards}      </div>\n'
                        f"    </div>"
                    )
            if not date_sections:
                continue
            month_sections.append(
                f'  <div class="month-section" data-month-id="{month_id}" id="{month_id}">\n'
                f'    <div class="month-heading"><span class="toggle-arrow">▼</span>{month_name} {year}</div>\n'
                f'    <div class="month-content">\n'
                + "\n".join(date_sections) + "\n"
                f'    </div>\n'
                f"  </div>"
            )
        if not month_sections:
            continue
        year_id = f"year-{scantype}-{year}"
        year_sections.append(
            f'<div class="year-section" id="{year_id}">\n'
            f'  <div class="year-heading"><span class="toggle-arrow">▼</span>{year}</div>\n'
            f'  <div class="year-content">\n'
            + "\n".join(month_sections) + "\n"
            f"  </div>\n"
            f"</div>"
        )

    content = "\n".join(year_sections) if year_sections else '<p class="empty">No images found.</p>'
    return (
        f'<div class="scan-panel" id="panel-{scantype}">\n'
        f"{content}\n"
        f"</div>"
    )


def build_sidebar_nav(year_month_dates):
    """Build the sidebar year/month navigation."""
    html = ""
    for year, months in year_month_dates.items():
        month_links = ""
        for month in months:
            month_name = MONTH_NAMES[int(month) - 1][:3]
            # href uses the first scan panel id; JS will rewrite anchors on tab switch
            month_links += (
                f'    <a class="month-link" href="#panel-__SCAN__-{year}-{month}">'
                f"{month_name}</a>\n"
            )
        html += (
            f'<div class="year-group">\n'
            f'  <button class="year-toggle"><span>{year}</span>'
            f'<span class="arrow">▶</span></button>\n'
            f'  <div class="month-list">\n{month_links}  </div>\n'
            f"</div>\n"
        )
    return html


def generate_html(root: Path, title: str, subtitle: str) -> str:
    data = discover_images(root)

    # Collect all scan types
    all_scantypes = []
    for st in SCAN_ORDER:
        if any(st in v for v in data.values()):
            all_scantypes.append(st)
    for date_scans in data.values():
        for st in date_scans:
            if st not in all_scantypes:
                all_scantypes.append(st)

    all_dates = sorted(data.keys(), reverse=True)

    if not all_scantypes:
        return "<html><body><p>No quicklook images found.</p></body></html>"

    year_month_dates = group_dates_by_year_month(all_dates)

    # Build per-scantype lookup: {scantype: {date: [paths]}}
    by_scan = defaultdict(lambda: defaultdict(list))
    for date_str, scans in data.items():
        for scantype, paths in scans.items():
            by_scan[scantype][date_str] = paths

    first_scan = all_scantypes[0]
    sidebar_nav = build_sidebar_nav(year_month_dates).replace("__SCAN__", first_scan)

    tab_buttons = ""
    scan_panels = ""
    for i, st in enumerate(all_scantypes):
        active = "active" if i == 0 else ""
        tab_buttons += (
            f'  <button class="tab-btn {active}" data-target="panel-{st}">'
            f"{st.upper()}</button>\n"
        )
        panel = build_scan_panel(st, by_scan, year_month_dates)
        if i == 0:
            panel = panel.replace('class="scan-panel"', 'class="scan-panel active"', 1)
        scan_panels += panel + "\n\n"

    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    return HTML_TEMPLATE.format(
        title=title,
        subtitle=subtitle,
        generated=generated,
        sidebar_nav=sidebar_nav.rstrip(),
        tab_buttons=tab_buttons.rstrip(),
        scan_panels=scan_panels.rstrip(),
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate a quicklooks index.html")
    parser.add_argument(
        "-o", "--outdir",
        default="/data/processing/kepler/reading-general/quicklooks",
        help="Root quicklooks directory (default: %(default)s)",
    )
    parser.add_argument(
        "--title",
        default="University of Reading Ka-band Cloud Radar (Kepler) — Quicklooks",
        help="Page title",
    )
    parser.add_argument(
        "--subtitle",
        default="University of Reading · STFC Chilbolton Observatory · Data processed with kepler-radar-utils v1.2.0",
        help="Page subtitle",
    )
    args = parser.parse_args()

    root = Path(args.outdir).resolve()
    if not root.is_dir():
        raise SystemExit(f"Directory not found: {root}")

    html = generate_html(root, args.title, args.subtitle)
    out = root / "index.html"
    out.write_text(html, encoding="utf-8")
    print(f"Written: {out}")


if __name__ == "__main__":
    main()
