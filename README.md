# CCkitViewer

CCkitViewer is a standalone desktop project, moved from CCkit's `visualize/`
directory to `D:\Github_repositories\CCkitViewer`. Run the commands below from
the CCkitViewer project root. Its source, feature backends, examples and build
configuration do not require a CCkit checkout. The application and packaged
executable are named **CCkitViewer**. The entry point is `viewer.py`, with
`run_viewer.bat` and `run_viewer_exe.bat` as the source and standalone launchers.
`AGENTS.md` documents the standalone project's development instructions;
all relative paths and commands in it use this repository root.

The desktop application is organized around registered visualization features.
The first connected feature is **Bicrystal Builder**, the grain-boundary builder backed by
`core/bicrystal/bicrystal_backend.py`. The feature selector and right-hand card are
designed so future CCkit scripts can add their own names and parameter cards;
the currently unimplemented entries remain visibly reserved instead of being
treated as working modelers. See [AGENTS.md](AGENTS.md) for the extension
contract.

The Windows desktop viewer supports crystal structures and trajectories,
Bicrystal Builder grain-boundary modeling, Surface Builder, Interface Builder, and Wyckoff crystal
generation. The left side renders structures with PyVista/VTK, while the
right side displays the selected feature's parameter card.

For Bicrystal Builder's build and scan modes, the right side starts with built-in
parameters based on `example/Cu/gbin.json`, initially in **Scan CSL candidates**
mode. Select a structure to run with
these defaults; importing `gbin.json` is optional and replaces the current
configuration. Valid build changes
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

The **Bicrystal Builder mode** selector changes between **Build boundary** and **Scan CSL
candidates**. Build exposes the rotation axis, angle, target Sigma, and
symmetric GB plane. Scan exposes the angle start, stop, increment, and a list
of successful results. Build-only live adjustments (`dp1`, `dp2`, CNID grid,
GB spacing, vacuum, orthogonality and mirror) and the model summary are hidden
in scan mode. The mode cards size themselves to the visible content.
CSL policy, rotation axis, maximum Sigma, tolerances, and both modes' values
are kept in the same JSON document; switching mode changes only `mode`, and
editing scan controls preserves the hidden build values. A scan result can be
copied into build mode, and successful rows can be exported as
`GB_scan_results.csv`.

The scan card shows a fixed-width table for quick reading. Deformation values
are shown with two decimal places when that is informative, and very small
values retain scientific notation. The exported file remains a standard CSV.

Numeric and option controls consume mouse-wheel events, so scrolling over a
value or a dropdown neither changes it nor puts it into an editing state.
Wrapped descriptions and summaries use Qt's width-dependent height calculation
during layout. Their cards allow vertical growth, and the scroll area updates
its content height. Label resize handlers do not force a minimum height after
the following rows have already been positioned. Shared label behavior is
implemented in `viewer_widgets.py`.

The Bicrystal Builder card also exposes InterfaceMaster's `mirror` option. It is written
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
├── bicrystal/
│   └── bicrystal_backend.py
├── surface/
│   ├── surface_backend.py
│   └── surface_panel.py
├── interface/
│   ├── interface_backend.py
│   ├── pymatgen_engine.py
│   ├── interface_master_engine.py
│   ├── commensurability.py
│   ├── registry_analysis.py
│   ├── interface_overlay.py
│   └── interface_panel.py
└── wyckoff/
    ├── wyckoff_model.py
    └── wyckoff_panel.py
```

The Bicrystal Builder backend was copied from CCkit's `Scripts/model/model_gb.py`;
its CLI original remains in CCkit. The GUI copy uses a distinct filename so
the CCkit command dispatcher's recursive search stays unambiguous.
Interface Builder dispatches to independent pymatgen and interface-master engines
under `core/interface/`. The InterfaceMaster adapter adds planar rational-period
analysis without changing the upstream checkout or using Bicrystal Builder's settings.
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
python -m pip install -r requirements-viewer.txt
```

InterfaceMaster is loaded from the local checkout in
`D:\Github_repositories\interface_master`; it is not installed from the
configured PyPI mirror. `run_viewer.bat` adds that sibling checkout to
`PYTHONPATH` automatically. If you run the Python file directly, set it first:

```powershell
$env:PYTHONPATH = "D:\Github_repositories\interface_master;$env:PYTHONPATH"
```

VTK uses the available OpenGL hardware renderer. A discrete GPU is helpful for
large supercells, but the viewer remains usable with software rendering.

## Start

From the repository root, pass a structure to use the built-in parameters:

```powershell
.\run_viewer.bat example\Cu\POSCAR_Cu
```

Optionally pass a configuration file to start with saved parameters:

```powershell
.\run_viewer.bat example\Cu\POSCAR_Cu example\Cu\gbin.json
```

Or run the Python entry point directly:

```powershell
python viewer.py example\Cu\POSCAR_Cu example\Cu\gbin.json
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

The Bicrystal Builder build card exposes **CNID u1 subdivisions** and **CNID u2
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

When no arguments are provided, the source viewer looks for `ZrO2.cif` in
`example/` and uses the built-in Bicrystal Builder parameters. JSON is loaded only when
passed as an argument or imported in the UI. Double-clicking the standalone exe opens
a viewer with default Bicrystal Builder parameters and file-loading buttons; choosing
a structure is sufficient to start without importing JSON. A brief splash
precedes the main window. Startup errors appear in a dialog, with a full traceback in
`%TEMP%\CCkitViewer_startup.log`. The viewer supports both `gbin.json`
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
powershell -NoProfile -ExecutionPolicy Bypass -File .\build_viewer.ps1 `
  -PythonExecutable D:\miniconda3\envs\interface\python.exe
.\run_viewer_exe.bat example\Cu\POSCAR_Cu example\Cu\gbin.json
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
`dist\CCkitViewer\CCkitViewer.exe` together with its runtime
directory. Keep the whole `CCkitViewer` directory when copying it to another
Windows computer.

Window titles, the startup splash, and error dialogs use the CCkitViewer name.
Existing display preferences and remembered folders are retained across the
rename. Rebuild with `build_viewer.ps1` to update an existing standalone
installation; previously built executables keep their original name and UI.


## Three-index and hexagonal four-index notation

Bicrystal Builder has independent notation selectors on the right of **Rotation axis**
and **Symmetric GB plane**. **Miller indices** use [UVW] for directions and (hkl)
for planes; **Miller–Bravais indices** use [uvtw] and (hkil). The redundant t/i index is
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

## Surface Builder: oriented slabs and atomic terminations

Select **Surface Builder** and import a bulk structure. Enter the surface plane
and two nonparallel in-plane crystal directions, all in the input cell's basis.
Three-index and standard hexagonal Miller–Bravais notation are supported.
The directions must satisfy the zone law `hU + kV + lW = 0` and form a
right-handed pair with the selected plane normal. The program reports invalid
directions instead of projecting or silently changing them.

Each direction is reduced to coprime integer indices in the input basis;
**Repeats along a/b** scale those translation vectors. The bulk is not
converted to another conventional or primitive cell. The generated surface
may therefore be larger than a minimal primitive surface cell. No strain,
surface reconstruction, passivation or atom deletion for close contacts is applied.
Output a/b follow the selected directions, with a along Cartesian x and b in
the xy plane. The c vector is perpendicular to the surface along z.

**Cut-window thickness** is the width of a bulk interval retained along the
normal. Since atomic layers are discrete, the outermost-atom thickness can
be smaller; the summary reports both. **Vacuum gap** is the exact normal
separation between outermost atomic planes of neighboring periodic slabs,
split equally above and below the centered slab.

Click **Generate surface** once to populate **Interlayer cut** choices. Each
choice cuts through an atomic-layer gap and lists its phase and bottom-layer
elements. Selecting another cut changes the exposed layers; the summary shows
the actual top and bottom elements. The cuts are not symmetry-filtered.
Alternatively, choose **Custom cut position** and specify the lower boundary
as a fraction of the reduced (hkl) normal repeat. The retained interval is
`[cut, cut + thickness)`; phase 1 is equivalent to phase 0 up to a bulk
translation. Top and bottom terminations are coupled to the chosen thickness.
**Layer tolerance** groups nearby atomic heights for cut enumeration and
exposed-element reporting; it does not move atoms.

The independent card retains its controls across feature switches. Optional
automatic rebuilding is initially off. The main viewer renders both surface
planes, preserves the camera on rebuild and exports POSCAR, CIF or XYZ.
Load/save **surface.json** through the card or the shared configuration actions;
relative bulk paths resolve beside the JSON file.

For example, load `example/Cu/POSCAR_Cu`, use plane `(001)` and directions
`[100]`, `[010]`, then set thickness and vacuum. A matching JSON is included
as `example/surface_Cu001.json`.

## Interface Builder: two-material heterostructures

Select **Interface Builder**, then choose the **Engine** directly below the
feature selector. **pymatgen** implements Zur and McGill's *Lattice match: an
application to heteroepitaxy* (1984). **interface-master** adds explicit rotation
searches, minimal planar CSL/DSCL analysis and rational-period diagnostics inspired
by Inoue et al. (2026), using InterfaceMaster's native lattice tools.

The following workflow describes **pymatgen**:

1. Select substrate and film files and their contact planes in each input
   cell's basis. Three indices and standard hexagonal four indices are supported.
2. Adjust maximum area, area-ratio filter, relative length/angle tolerances,
   bidirectional comparison, candidate limit, and termination enumeration.
3. Click **Search matching cells**. Searches run in a background thread and
   can be cancelled between native search steps. Retained results are sorted
   by film area, then strain. A truncated candidate list is explicitly marked.
4. Select a matching cell. Compare both unstrained areas, signed length and
   angle mismatches, film von Mises strain, actual supercell matrices, and
   reduced lattice vectors. The diagram overlays the two reduced cells at
   equal scale with their first vectors aligned; it is not an atomic preview
   or a crystallographic twist-angle plot.
5. Choose a **film / substrate termination** pair. Set each slab's minimum
   thickness in oriented unit planes or Å, the interface gap and external vacuum.
6. Choose **Custom shift** for one structure at specified fractional a/b
   offsets, or **Grid sampling** for multiple shifts starting at (0, 0).
   Only the active mode's controls are shown; switching preserves both sets
   of values. Click **Build selected interface**. View the atomic model in the main 3D viewport and
   export the current structure as POSCAR/CIF/XYZ or all shift frames as XYZ.

The implementation uses pymatgen's `ZSLGenerator` and
`CoherentInterfaceBuilder.get_interfaces`. The substrate is the reference;
the film's in-plane lattice is strained into coherency. Shift sampling covers
the matched supercell's a/b vectors, not a CNID cell. Geometry edits reuse
cached lattice searches; changing files, planes, or matching settings requires
a new search. Auto rebuilding selected geometry is optional and initially off.
All numeric and selection controls consume mouse-wheel events.

**Load config…** routes to this card. The collapsible **interface.json** group
also provides load/save/apply actions. Relative structure paths in loaded JSON
resolve beside the JSON file; a manually entered relative film path resolves
beside the substrate. Old InterfaceMaster JSON settings are rejected explicitly
instead of silently acquiring different semantics.

For **interface-master**, choose the rotation interval/step, rational denominator,
principal strain, area-change, CSL-index and commensurability limits, then search.
The finite mesh is scanned completely; the retained list is ranked by σ, strain
or area. A second diagram shows ideal planar lattice sites and a DSCL cell;
continued-fraction and Farey guides explain 1D matching periods. Select the
terminations, thickness, gap and vacuum, then build the atomic structure.

Both interface engines keep the two in-plane vectors as lattice a/b and the
stacking/normal direction as c in previews and exports. The common plane
renderer reads `interface_normal_axis=2` from frame metadata, so interface
markers and the DSCL overlay follow fractional c without reordering the cell.
Older frames without this marker retain their first-axis convention; Bicrystal
Builder keeps its existing convention. This metadata stays internal and does
not add custom properties to the XYZ comment line.

For **interface-master**, custom shifts can use either **DSCL cell** vectors or
the **Matched a/b supercell** vectors. Enter coefficients **d1** and **d2** to
generate one structure at `d1 × vector1 + d2 × vector2`. The new setting
`custom_shift_basis` defaults to `supercell`, preserving old custom-shift JSON
semantics. Grid sampling retains its independent `im_shift_basis` selection.
The viewport status line shows `d1` and `d2` as coefficients of the selected
shift-cell vectors and names the basis (DSCL cell or matched supercell).
These are the original coefficients, before conversion to matched-supercell
fractional coordinates or periodic wrapping. For a 2×3 grid, the last frame
(index 5) therefore shows `d1=0.5, d2=0.7`, while calculations still use
`1/2` and `2/3`. Older trajectories containing
only `in_plane_offset` label those values as supercell fractional shifts.
For **pymatgen**, coefficients always refer to the matched a/b supercell.
Each built frame reports
projected atomic registry, nearest cross-interface distance and close contacts;
the frame slider updates these diagnostics. They ignore species in the distance
metric and do not represent interface energies. No atom is silently removed.
Built InterfaceMaster models also show a purple DSCL outline and its two
translation vectors on the interface plane in the left 3D view, at their actual
lengths. This overlay follows the plane visibility preference and the current
frame. The right-side diagram retains a separately scaled DSCL inset.
**Save matching / model report…** writes search results, matrices, matching-period
guides and available current-model/frame diagnostics to JSON.

Enter a **System label** (for example `hZrO-Zr`) in the atomic-interface card;
an empty label uses the substrate-film reduced formulas. Export **Current
frame** or **All frames** as **XYZ · extended format (*.xyz)** to retain the
cell, PBC and per-frame `config_type` in the extended comment line. For
interface-master the format is:

```text
hZrO-Zr_sigma_6.00_theta_30.0_zspacing_2.0_d1_0.5_d2_0.3
```

Sigma and theta refer to the actual candidate, formatted with two and one
decimal places respectively in `config_type` and the interface result displays;
calculations retain the original precision. zspacing is the interface gap
(`gap`) in Å, and d1/d2 are coefficients in the selected basis. These three
fields use exactly one decimal in `config_type`, their numeric controls and
frame summaries. For example, a coefficient of `0.333333333333` is displayed
as `0.3`; its full value remains in calculations and JSON reports. Editing or
collecting an unchanged control does not replace its value with the rounded
display. This rule applies to both engines. Interface XYZ
exports contain only `config_type` as a custom property, alongside the standard
`Lattice`, `Properties` and `pbc` fields. Shift vectors, frame indices,
terminations, registry and other diagnostic metadata remain in the viewer and
JSON report; they are not written to the XYZ comment line.
Pymatgen uses `{system}_zsl_match_{number}_strain_{value}_zspacing_{value}_d1_{value}_d2_{value}`
with the candidate's dimensionless von Mises strain, without inventing a Sigma
or absolute twist angle. The current frame's tag is displayed for copying.

See [pymatgen parameters and methodology](docs/interface_pymatgen.md) and
[InterfaceMaster/Inoue methodology](docs/interface_master.md) for definitions,
defaults and limitations. Source and PyInstaller configuration include both
engines. Rebuild manually to update the standalone program.

## Wyckoff model: build crystals from symmetry positions

Select **Wyckoff model** from **Visualization feature**. The independent card
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

PyXtal is included in `requirements-viewer.txt`; install these requirements
in the Python environment used to launch or build the viewer. The PyInstaller
configuration includes both Wyckoff modules and PyXtal's database files.
Source changes require a manual rebuild to appear in the standalone `.exe`.
