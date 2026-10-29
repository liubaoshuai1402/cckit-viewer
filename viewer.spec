# PyInstaller onedir build for the standalone CCkitViewer application.
from pathlib import Path

import six
from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_submodules


VISUALIZE_DIR = Path(SPECPATH)
REPOSITORY_ROOT = VISUALIZE_DIR
INTERFACE_MASTER_ROOT = REPOSITORY_ROOT.parent / "interface_master"

datas = [
    # six is a single-file module.  The current PyInstaller hooks can analyze
    # its dynamic ``six.moves`` imports without copying the source file, but
    # the frozen PyVista stack still needs the module at runtime.
    (str(Path(six.__file__)), "."),
]
# Backends are loaded by file path. Retain the same feature-directory layout
# below _MEIPASS as in the source tree, including future core modules.
for core_source in sorted((VISUALIZE_DIR / "core").rglob("*.py")):
    destination = core_source.parent.relative_to(VISUALIZE_DIR).as_posix()
    datas.append((str(core_source), destination))
for sample_name in ("ZrO2.cif", "gbin.json", "surface_Cu001.json", "Cu/POSCAR_Cu"):
    sample_path = REPOSITORY_ROOT / "example" / sample_name
    if sample_path.is_file():
        datas.append((str(sample_path), sample_path.parent.relative_to(REPOSITORY_ROOT).as_posix()))

binaries = []
hiddenimports = [
    "viewer_widgets",
    "six",
    "ase",
    "ase.io.extxyz",
    "ase.io.trajectory",
    "pyxtal.symmetry",
    "core.bicrystal.bicrystal_backend",
    "core.surface.surface_backend",
    "core.surface.surface_panel",
    "core.interface.interface_backend",
    "core.interface.interface_panel",
    "core.interface.interface_overlay",
    "core.interface.pymatgen_engine",
    "core.interface.interface_master_engine",
    "core.interface.commensurability",
    "core.interface.registry_analysis",
    "scipy.spatial",
    "pymatgen.analysis.interfaces.zsl",
    "pymatgen.analysis.interfaces.coherent_interfaces",
    "pymatgen.analysis.elasticity.strain",
    "pymatgen.core.interface",
    "pymatgen.core.surface",
    "core.wyckoff.wyckoff_model",
    "core.wyckoff.wyckoff_panel",
    "interfacemaster",
    "interfacemaster.cellcalc",
    "interfacemaster.interface_generator",
    "interfacemaster.tool",
]
# PyInstaller ships maintained hooks for NumPy, PySide6, PyVista, and VTK.
# Collect only the small packages whose imports are dynamic in this project;
# collecting every VTK and Qt submodule makes the build needlessly slow.
for package_name in ("pyvistaqt", "interfacemaster"):
    package_datas, package_binaries, package_hiddenimports = collect_all(package_name)
    datas.extend(package_datas)
    binaries.extend(package_binaries)
    hiddenimports.extend(package_hiddenimports)

# pymatgen loads its element table and related compressed data at runtime.
datas.extend(collect_data_files("pymatgen"))
datas.extend(collect_data_files("pyxtal"))

hiddenimports.extend(collect_submodules("interfacemaster"))

# The build host has another conda environment on PATH.  Its ICU DLLs are
# incompatible with the Qt ICU imports and would shadow Windows' own ICU at
# runtime.  Qt can use the system ICU on supported Windows versions.
def _is_icu_file(item: tuple[str, str]) -> bool:
    """Return whether a collected file belongs to the ICU runtime."""

    return Path(item[0]).name.lower().startswith("icu")


datas = [item for item in datas if not _is_icu_file(item)]
binaries = [
    item for item in binaries if not _is_icu_file(item)
]

analysis = Analysis(
    [str(VISUALIZE_DIR / "viewer.py")],
    pathex=[str(VISUALIZE_DIR), str(REPOSITORY_ROOT), str(INTERFACE_MASTER_ROOT)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)
pyz = PYZ(analysis.pure)
exe = EXE(
    pyz,
    analysis.scripts,
    [],
    exclude_binaries=True,
    name="CCkitViewer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
)


def _is_runtime_icu(entry: tuple[str, ...]) -> bool:
    """Exclude ICU files added by PyInstaller's automatic Qt hooks."""

    return any(Path(str(value)).name.lower().startswith("icu") for value in entry[:2])


runtime_binaries = [
    entry for entry in analysis.binaries if not _is_runtime_icu(entry)
]
runtime_datas = [entry for entry in analysis.datas if not _is_runtime_icu(entry)]

COLLECT(
    exe,
    runtime_binaries,
    runtime_datas,
    strip=False,
    upx=True,
    name="CCkitViewer",
)
