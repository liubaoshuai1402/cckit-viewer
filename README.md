# CCkitViewer

CCkitViewer is a standalone desktop project, moved from CCkit's `visualize/`
directory to `D:\Github_repositories\CCkitViewer`. Run the commands below from
the CCkitViewer project root. Its source, feature backends, examples and build
configuration do not require a CCkit checkout. The existing executable and
launcher filenames are retained. `AGENTS.md` is preserved unchanged from the
original directory; its `visualize/` paths refer to this project root, and its
references to CCkit documentation describe the original repository layout.

The desktop application is organized around registered visualization features.
The first connected feature is **OmniGB**, the grain-boundary builder backed by
`core/omnigb/omnigb_backend.py`. The feature selector and right-hand card are
designed so future CCkit scripts can add their own names and parameter cards;
the currently unimplemented entries remain visibly reserved instead of being
treated as working modelers. See [AGENTS.md](AGENTS.md) for the extension
contract.

This folder contains a Windows desktop viewer for the OmniGB `model-gb.py`
build and scan modes. The left side renders a generated bicrystal with
PyVista/VTK. The right side loads and edits `gbin.json`; valid build changes
rebuild the preview after a short debounce, while scan mode runs the CSL angle
sweep in a worker and displays the successful rows. The `dp1` and `dp2`
termination controls have their own adjustment steps; the two CNID subdivision
controls determine the RBT frame grid.
The top **Personal preferences** toolbar opens settings for atom size, bond
visibility and cutoffs, and symmetric GB plane visibility, opacity, and color.
The two displayed planes are the two periodic boundaries of the final
bicrystal cell. The viewport has a white background with alpha transparency.
Both the orientation marker and the plot axes explicitly label x, y, and z;
the labels are enlarged and bold, with red, green, and blue direction arrows.

The **OmniGB mode** selector changes between **Build boundary** and **Scan CSL
candidates**. Build exposes the rotation axis, angle, target Sigma, and
symmetric GB plane. Scan exposes the angle start, stop, increment, and a list
of successful results. CSL policy, rotation axis, maximum Sigma, tolerances,
and all live adjustment values are kept in the same JSON document; switching
mode changes only `mode`, so shared values are not reset. A scan result can be
copied into build mode, and successful rows can be exported as
`GB_scan_results.csv`.

The scan card shows a fixed-width table for quick reading. Deformation values
are shown with two decimal places when that is informative, and very small
values retain scientific notation. The exported file remains a standard CSV.

Numeric and option controls consume mouse-wheel events, so scrolling over a
value or a dropdown neither changes it nor puts it into an editing state.

The OmniGB card also exposes InterfaceMaster's `mirror` option. It is written
as the JSON boolean `mirror` in `gbin.json`; enabling it asks InterfaceMaster
to construct the second grain from a reflected copy of the first grain.
The same card exposes `GB spacing` (`normal_distance`) and `External vacuum`
(`vacuum`), both in Angstrom. They change the cell geometry without adding
atoms. The `normal_ortho` and `plane_ortho` checkboxes control whether
InterfaceMaster searches for approximately normal and approximately
orthogonal CSL cell vectors.

## Feature implementation directories

Feature code is organized under `core/`:

```text
core/
├── omnigb/
│   └── omnigb_backend.py
├── interface/
│   ├── interface_backend.py
│   └── interface_panel.py
└── wyckoff/
    ├── wyckoff_model.py
    └── wyckoff_panel.py
```

The OmniGB and Interface Builder backends were copied from
CCkit's `Scripts/model/model_gb.py` and `Scripts/model/model_heterointerface.py`.
Their CLI originals remain in CCkit. The GUI copies use distinct
filenames so the CCkit command dispatcher's recursive search stays unambiguous.
Interface Builder locates shared OmniGB helpers in `core/omnigb/`.
The Wyckoff backend and the independent parameter card modules have moved into
their feature directories. The application shell and registry remain at the
top level. New feature implementations must follow this organization.

The registry loads the core backends for both source and frozen builds.
PyInstaller includes their source files under `core/` in its runtime
directory and collects the feature packages and their dependencies. Existing
sample inputs are read from this project's `example/` directory.

## Install

Use a build environment that contains pymatgen, the GUI packages, and
PyInstaller:

```powershell
python -m pip install -r requirements-gb-viewer.txt
```

InterfaceMaster is loaded from the local checkout in
`D:\Github_repositories\interface_master`; it is not installed from the
configured PyPI mirror. `run_gb_viewer.bat` adds that sibling checkout to
`PYTHONPATH` automatically. If you run the Python file directly, set it first:

```powershell
$env:PYTHONPATH = "D:\Github_repositories\interface_master;$env:PYTHONPATH"
```

VTK uses the available OpenGL hardware renderer. A discrete GPU is helpful for
large supercells, but the viewer remains usable with software rendering.

## Start

From the repository root, pass a structure and configuration file:

```powershell
.\run_gb_viewer.bat example\Cu\POSCAR_Cu example\Cu\gbin.json
```

Or run the Python entry point directly:

```powershell
python gb_viewer.py example\Cu\POSCAR_Cu example\Cu\gbin.json
```

The top toolbar also has **Load structure…** and **Load config…** actions,
so the standalone program can be opened first and files selected afterward.
Open dialogs and CSV/JSON save dialogs remember the last selected folder
across application restarts. Structure and trajectory exports default to the
source file's directory.
The structure picker accepts CIF, POSCAR/CONTCAR, XYZ, PDB, RES, CSSR, XSF,
MSON/JSON, and YAML files; **All files** remains available for other formats
understood by pymatgen. Use **Load trajectory…** for multi-frame extended XYZ
(`.xyz` or `.extxyz`) and ASE `.traj` files. Every lattice trajectory frame must
contain a valid cell; coordinate-only XYZ files do not contain enough lattice
information for the viewer.

Use **Personal preferences** to set the global bond range. After a preview
has been built, the same dialog lists every element pair in the structure;
each pair can be enabled or disabled and given its own minimum and maximum
distance. Pairs not listed in a saved preference use the global range.

## CNID grid and frame navigation

The OmniGB build card exposes **CNID u1 subdivisions** and **CNID u2
subdivisions** in place of the Cartesian translation controls. The JSON key
is `"cnid_grid": [n1, n2]`; the viewer defaults to `[4, 4]`. Each direction
accepts a positive integer up to 1000. The grid generates `n1*n2` frames with

```text
translation = (i/n1)*u1 + (j/n2)*u2
frame = i*n2 + j
i = 0,...,n1-1; j = 0,...,n2-1
```

The endpoint at coefficient 1 is excluded because it repeats the periodic
boundary. CSL and CNID are computed once for the grid; every frame keeps the
same `dp1`, `dp2`, normal spacing, vacuum and supercell settings. The legacy
`translation` and `dydz` keys are removed from GUI configurations. The public
`CCkit model-gb` command retains its existing single-POSCAR behavior and legacy
Cartesian translation support; `cnid_grid` controls the GUI sampling sequence.
For `mirror: true`, the CNID basis uses the actual mirrored grains' common
surface translation lattice.

The viewport has a horizontal frame slider and a synchronized integer editor
to its right. Both use zero-based indices: the first frame is 0 and the last
is `n1*n2-1`. Changing frames preserves the camera and does not rebuild the
model. The status line shows the selected frame and its CNID coefficients.
Grid generation runs in a worker thread and reports progress; superseded
previews stop between frames so the newest settings can be built.

## Export a frame or the complete sequence

**Export structure / trajectory…** opens a dialog to select **Current frame**
or **All frames**. A current frame can be written as POSCAR, CIF or extended
XYZ. All frames are written as one multi-frame extended XYZ file. The `.xyz`
output retains each frame's lattice, PBC, frame index, CNID coefficients,
grid indices and translation vector. It can be reopened with **Load
trajectory…**. These are independent RBT candidates in a deterministic order.

Exports default to the source structure or trajectory directory. The export
captures the settings, frame index and sequence when clicked; later frame
changes do not change the saved selection. If current model parameters have
not been built yet, the matching grid is built in the background. A failed
build never exports the older preview. Writing also runs in the background;
the destination is replaced only after the complete file has been written.

Changing `gbin.json` keeps the current VTK camera position, direction, and
zoom. The first preview establishes the initial view; later rebuilds reuse it.

When no arguments are provided, the source viewer looks for `ZrO2.cif` and
an optional `gbin.json` in `example/`. Double-clicking the standalone exe opens
an empty viewer with file-loading buttons. A brief splash precedes the main
window. Startup errors appear in a dialog, with a full traceback in
`%TEMP%\CCkitGBViewer_startup.log`. The viewer supports both `gbin.json`
`mode: "build"` and `mode: "scan"`.

The preview is built in a temporary directory and does not overwrite the
POSCAR or grain CIF files beside the input. Use **Save gbin.json** when the
edited values should be written to the loaded configuration, or **Save current
JSON as…** to write the live editor state to a chosen file. The default name
for JSON saves is `gbin.json`.

## Build a standalone Windows program

The build uses PyInstaller's onedir layout. It bundles Python, Qt, VTK,
PyVista, and the local InterfaceMaster code; the resulting executable does
not require Python or InterfaceMaster on the target machine.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\build_gb_viewer.ps1 `
  -PythonExecutable D:\miniconda3\envs\interface\python.exe
.\run_gb_viewer_exe.bat example\Cu\POSCAR_Cu example\Cu\gbin.json
```

If the `interface` environment is already active, `-PythonExecutable` can be
omitted.

Keep `build` to reuse PyInstaller's build cache. Add `-Clean` only
when a dependency change requires a full rebuild. A build must be able to read
all files in the build environment; a restricted build can omit required
modules even when PyInstaller exits successfully.

The spec bundles `six.py` and pymatgen's runtime data explicitly. It also
excludes conda ICU DLLs that would shadow Windows' system ICU and prevent
QtCore from loading. Validate the frozen app itself with an exposed native
window and an actual preview build, rather than checking only process lifetime.

The portable program is created at
`dist\CCkitGBViewer\CCkitGBViewer.exe` together with its runtime
directory. Keep the whole `CCkitGBViewer` directory when copying it to another
Windows computer.


## Three-index and hexagonal four-index notation

OmniGB has independent notation selectors on the right of **Rotation axis**
and **Symmetric GB plane**. **三指数** uses [UVW] for directions and (hkl)
for planes; **六方四指数** uses [uvtw] and (hkil). The redundant t/i index is
computed automatically and displayed read-only. Changing the selector
converts the values rather than reinterpreting the same digits.

JSON settings use `rotation_axis_notation` and `gb_plane_notation`, each set
to `"three"` or `"four"`. Existing three-index files remain accepted; absent
notation keys are inferred from vector length. For example, the hexagonal
c axis and basal plane are:

```json
{
  "rotation_axis": [0, 0, 0, 1],
  "rotation_axis_notation": "four",
  "gb_plane": [0, 0, 0, 1],
  "gb_plane_notation": "four"
}
```

Four-index input requires a standard hexagonal input cell: a=b,
gamma=120 degrees, and c perpendicular to the basal plane. A rhombohedral
primitive setting or a rotated/redefined nonstandard lattice basis must use
three indices in its actual basis or be converted to a standard hexagonal
cell first. The program does not silently standardize input structures.
Directions convert to [u-t, v-t, w]; planes convert to (h,k,l). Conventional
to primitive-cell conversion is then performed by InterfaceMaster.

## Interface Builder: two-material heterostructures

Select **Interface Builder** from **Visualization feature**. Its independent
card uses `core/interface/interface_backend.py` and owns `interface.json`;
OmniGB's input and JSON are retained when switching cards. The existing
`model_interface.py` coherent-interface command is unaffected.

1. Choose a substrate and a film structure.
2. Specify each contact plane in that material's own input-cell basis. Each
   plane has its own three/four-index selector, permitting cubic/hexagonal pairs.
3. Optionally enable reference-direction alignment and enter one in-plane
   direction for each material, again with separate notation selectors.
4. Set the twist-angle interval and step, maximum substrate/film surface
   area indices, and matching/deformation tolerances.
5. Set separate normal repeats, terminations, spacing, vacuum and CNID grid.
6. Choose **Find match and build CNID frames**. Automatic updates can be disabled.

The backend uses native `search_one_position_2D`, `compute_bicrystal_two_D`,
and `get_bicrystal(two_D=True)` calls. It takes the first admissible 2D CSL
match in ascending twist order, retaining the substrate lattice and applying
the accepted transform to the film. It does not rank all possible matches
by energy. Equal start/stop angles select one twist angle. Search stops
between angles or sampling frames when newer settings supersede the request.

The result summary reports twist, both surface area indices, area, actual
film metric strain and atom count. Exact mode follows InterfaceMaster's
native exact-CSL policy; it must not be interpreted as a guarantee of an
unchanged film metric. Near mode makes the accepted deformation explicit.

The same zero-based frame slider and current/all-frame exports are shared
with OmniGB. Interface plane markers follow the separate slab thicknesses
and spacing rather than assuming the junction lies halfway through the cell.
With external vacuum, only the internal junction is marked.

**Load config…** routes to the selected card. The Interface Builder card also
provides **Load JSON…** and **Save JSON…**. Relative structure paths in a loaded
JSON file resolve beside that JSON; a manually entered relative film path
resolves beside the substrate. Export defaults to the substrate directory.

Both the new backend and the independent card are included in the PyInstaller
configuration. Rebuild the application manually to include these source changes
in the standalone program.

## wyckoff model: build crystals from symmetry positions

Select **wyckoff model** from **Visualization feature**. The independent card
uses `core/wyckoff/wyckoff_model.py` and does not require an input structure file:

1. Enter a space-group number from 1 to 230. Its Hermann–Mauguin symbol updates
   automatically, together with the available Wyckoff positions.
2. Enter the independent lattice constants in Å. Cubic groups show `a`;
   tetragonal, trigonal and hexagonal groups show `a, c`; orthorhombic groups
   show `a, b, c`; monoclinic groups additionally show `beta`; triclinic groups
   show all three lengths and all three angles in degrees. Each lattice
   parameter is displayed on its own row.
3. Use **Add Wyckoff position** to add or remove sites. Each numbered site card
   has its own border, an element symbol and a Wyckoff label such as `2c`.
   Selecting a position computes all equivalent coordinates in one orbit and
   displays them in **Select equivalent coordinates**. Choose one as the
   representative; generating the crystal expands the complete orbit, so
   choosing another equivalent representative gives the same crystal.
   Fixed positions list their exact coordinates. Positions with free parameters
   show all symbolic expressions with numerical values for one example orbit;
   this finite list cannot enumerate the infinitely many parameter values.
   The card explicitly identifies this example and the number of free parameters.
   For parameterized positions, **Free parameters (fractional)** appears below
   the equivalent-coordinate selector, with a separate row for each variable
   used by the symbolic expressions. For example, `(x, 0, 0)` shows only `x`.
   Enter a decimal or fraction and press Enter or leave the field to update all
   equivalent coordinates while retaining the selected expression. Generating
   or saving also commits these parameter edits. Invalid values are displayed
   inside the card and cannot generate or save an older orbit.
   Choose **Manual input** to reveal the fractional `x, y, z` editors. These
   accept decimals and simple fractions such as `1/3`. Switching back to
   **Select equivalent coordinates** computes the orbit from the supplied
   coordinates; invalid inputs keep the manual editors open with an error.
   The coordinate constraint and multiplicity remain visible in both modes.
   Changing the space group resets each site's position and coordinates to the
   new group's general position while retaining its element.
4. Click **Generate crystal** to expand the equivalent atoms and show the
   conventional cell in the left viewport. Edits require another click to
   update the preview. Coordinates inconsistent with the selected position,
   special coordinates with a different multiplicity, invalid cells and
   overlapping orbits produce an error instead of a partial structure.

The default example is hcp Zr: space group 194 (`P6_3/mmc`), `a=3.231`,
`c=5.147`, Zr on `2c` at `(1/3, 2/3, 1/4)`, producing two atoms.
The implementation follows the Wyckoff expansion used in
`undeveloped/hcp_zr_wyckoff.py`, with PyXtal supplying the symmetry tables and
pymatgen supplying the final structure. It uses spglib's default conventional
settings consistently for the displayed coordinates and symmetry operations:
monoclinic unique axis **b**, and trigonal **R** groups in hexagonal axes.
Alternative origin and axis settings are not selectable. Sites are fully
occupied; coincident species and partial occupancies are not supported.

**Load JSON…** and **Save JSON…** preserve the space group, lattice and site
list and each site's `coordinate_mode` in `wyckoff.json`; the toolbar's
**Load config…** action uses the selected card. Existing JSON files without
`coordinate_mode` restore their coordinates in manual mode.
After generating the current parameters, **Export structure / trajectory…**
can save POSCAR, CIF or extended XYZ. With no input file, the export dialog
uses the remembered folder (or the loaded configuration's directory as a
fallback). Generated bulk crystals do not display grain-boundary planes.

PyXtal is included in `requirements-gb-viewer.txt`; install these requirements
in the Python environment used to launch or build the viewer. The PyInstaller
configuration includes both Wyckoff modules and PyXtal's database files.
Source changes require a manual rebuild to appear in the standalone `.exe`.
