#! /usr/bin/env python

import os
import sys
import pathlib
import inspect

import waves
import yaml
from SCons.Node.Alias import default_ans

# Accept command line options with fall back default values
AddOption(
    "--build-dir",
    dest="variant_dir_base",
    default="build",
    nargs=1,
    type="string",
    action="store",
    metavar="DIR",
    help="SCons build (variant) root directory. Relative or absolute path. (default: '%default')"
)
AddOption(
    "--ignore-documentation",
    dest="ignore_documentation",
    default=False,
    action="store_true",
    help="Boolean to ignore the documentation build, e.g. during Conda package build and testing. Unaffected by the " \
         "'--unconditional-build' option. (default: '%default')"
)
AddOption(
    "--solve-cpus",
    dest="solve_cpus",
    default=1,
    nargs=1,
    type="int",
    action="store",
    metavar="N",
    help="Run the Abaqus Solve task using N CPUs. (default: '%default')"
)
AddOption(
    "--filter",
    dest="filter",
    default=False,
    action="store_true",
    help="Boolean speciyfing whether or not to run micromorphic filter for a particular upscaling study."\
         " (default: '%default')"
)
AddOption(
    "--calibrate",
    dest="calibrate",
    default=False,
    action="store_true",
    help="Boolean speciyfing whether or not to run calibration for a particular upscaling study."\
         " (default: '%default')"
)
AddOption(
    "--macro",
    dest="macro",
    default=False,
    action="store_true",
    help="Boolean speciyfing whether or not to run macro simulation in Tardigrade-MOOSE."\
         " (default: '%default')"
)
AddOption(
    "--macro-plastic",
    dest="macro_plastic",
    default=False,
    action="store_true",
    help="Boolean speciyfing whether or not to run heterogeneous, plastic macro simulation in Tardigrade-MOOSE. TEMPORARY! TO be replaced with dedicated Workflow in the future"\
         " (default: '%default')"
)
AddOption(
    "--macro-ignore-BCs",
    dest="macro_ignore_BCs",
    default=False,
    action="store_true",
    help="Boolean speciyfing whether or not to run macro simulation in Tardigrade-MOOSE"\
         " with an 'averaged' material card for boundary elements(default: '%default')"
)
AddOption(
    "--summary",
    dest="summary",
    default=False,
    action="store_true",
    help="Boolean speciyfing whether or not to run summary post-processing for multi-domain studies."\
         " (default: '%default')"
)
AddOption(
    "--peta-data-copy",
    dest="peta_data_copy",
    default=False,
    action="store_true",
    help="Boolean to create/update a local copy of the DNS data stored in the Peta Library. (default: '%default')"
)
AddOption(
    "--config-software",
    dest="config_software",
    default=False,
    action="store_true",
    help="Boolean to configure software paths. (default: '%default')"
)
# Inherit user's full environment and set project options
env = Environment(ENV=os.environ.copy(),
                  variant_dir_base=GetOption("variant_dir_base"),
                  ignore_documentation=GetOption("ignore_documentation"),
                  solve_cpus=GetOption("solve_cpus"),
                  filter=GetOption("filter"),
                  calibrate=GetOption("calibrate"),
                  macro=GetOption("macro"),
                  macro_plastic=GetOption("macro_plastic"),
                  macro_ignore_BCs=GetOption("macro_ignore_BCs"),
                  summary=GetOption("summary"),
                  peta_data_copy=GetOption("peta_data_copy"),
                  config_software=GetOption("config_software"),
                  TARFLAGS="-c -j",
                  TARSUFFIX=".tar.bz2"
)

# Empty defaults list to avoid building all simulation targets by default
env.Default()

# ============================================================ LINK SOFTWARE ===
# Sbatch
env['sbatch'] = waves.scons_extensions.add_program(["sbatch"], env)

# Sphinx
env["sphinx_build"] = waves.scons_extensions.add_program(["sphinx-build"], env)

# Read in config.yml
config_file = 'config_software.yml'
stream = open(config_file, 'r')
program_paths = yaml.load(stream, Loader=yaml.FullLoader)
stream.close()

# MPI
mpi_location = program_paths['mpi']
env['mpi'] = waves.scons_extensions.find_program(mpi_location, env)

# Abaqus
abaqus_windows = ["C:/Simulia/Commands/abaqus.bat"]
env["abaqus"] = waves.scons_extensions.find_program(program_paths['Abaqus'] + abaqus_windows, env)

# Cubit
env["cubit"] = waves.scons_extensions.add_cubit(program_paths['Cubit'] + ["cubit"], env)

# Ratel
ratel_location = program_paths['Ratel']
env['Ratel'] = waves.scons_extensions.find_program(ratel_location, env)
if env['Ratel']:
    env.PrependENVPath("PATH", str(pathlib.Path(env['Ratel']).parent))
ratel_solver = Builder(
    action=["cd ${TARGET.dir.abspath} && ${Ratel_program} \
             -options_file ${options_file} \
             -dm_plex_filename ${mesh_file} \
             -ts_monitor_diagnostic_quantities vtk:${monitor_file} \
             -ts_monitor_surface_force ascii:${force_file}:ascii_csv \
             -diagnostic_order 1 \
             > ${stdout_file} "])
ratel_solver_mpi = Builder(
    action=["cd ${TARGET.dir.abspath} && ${mpi_location} \
             -n ${ratel_cpus} ${Ratel_program} \
             -options_file ${options_file} \
             -dm_plex_filename ${mesh_file} \
             -ts_monitor_diagnostic_quantities vtk:${monitor_file} \
             -ts_monitor_surface_force ascii:${force_file}:ascii_csv \
             -diagnostic_order 1 \
             > ${stdout_file} "])
ratel_solver_sbatch = Builder(
    action=["cd ${TARGET.dir.abspath} && srun \
             -n ${ratel_cpus} ${Ratel_program} \
             -options_file ${options_file} \
             -dm_plex_filename ${mesh_file} \
             -ts_monitor_diagnostic_quantities vtk:${monitor_file} \
             -ts_monitor_surface_force ascii:${force_file}:ascii_csv \
             -diagnostic_order 1 \
             > ${stdout_file} "])
def ratel_builder_select():
    if env['sbatch']:
        return ratel_solver_sbatch
    elif env['mpi']:
        return ratel_solver_mpi
    else:
        return ratel_solver

# GEOS
## TODO: add GEOS

# Tardigrade Micromorphic Filter
filter_location = program_paths['filter'][-1]

# Calibration-micromorphic
micromorphic_location = program_paths['micromorphic'][-1]

# Calibration-constraints
constraints_location = program_paths['constraints'][-1]

# Tardigrade-MOOSE
tardigrade_location = program_paths['Tardigrade']
env['Tardigrade'] = waves.scons_extensions.find_program(tardigrade_location, env)
if env['Tardigrade']:
    env.PrependENVPath("PATH", str(pathlib.Path(env['Tardigrade']).parent))
tardigrade_solver = Builder(
    action=["cd ${TARGET.dir.abspath} && ${tardigrade_program} \
            -i ${tardigrade_input} \
            --n-threads=${tardigrade_cpus} \
            --no-color --color off > ${stdout_file} "])
tardigrade_solver_mpi = Builder(
    action=["cd ${TARGET.dir.abspath} && ${mpi_location} \
            -n ${tardigrade_cpus} ${tardigrade_program} \
            -i ${tardigrade_input} \
            --n-threads=2 \
            --no-color --color off > ${stdout_file} "])
tardigrade_solver_sbatch = Builder(
    action=["cd ${TARGET.dir.abspath} && srun \
            -n ${tardigrade_cpus} ${tardigrade_program} \
            -i ${tardigrade_input} \
            --n-threads=2 \
            --no-color --color off > ${stdout_file} "])

def tardigrade_builder_select():
    if env['sbatch']:
        return tardigrade_solver_sbatch
    elif env['mpi']:
        return tardigrade_solver_mpi
    else:
        return tardigrade_solver

# # Custom Paraview image generator
# env['paraview'] = waves.scons_extensions.find_program(program_paths['filter'], env)
# if env['paraview']:
    # env.PrependENVPath("PATH", str(pathlib.Path(env['paraview']).parent))
# paraview_image = Builder(
    # action=["cd ${TARGET.dir.abspath} && \
             # export MESA_GL_VERSION_OVERRIDE=3.2 && \
             # paraview --disable-xdisplay-test --force-offscreen-rendering --script ${script} ${script_options}"])

# ==================================================================== SCONS ===
# Set project internal variables and variable substitution dictionaries
project_configuration = pathlib.Path(inspect.getfile(lambda: None))
project_dir = project_configuration.parent
project_name = project_dir.name
version = "0.1.0"
author_list = ["Thomas Allard"]
author_latex = r" \and ".join(author_list)
latex_project_name = project_name.replace("_", "-")
documentation_source_dir = pathlib.Path("docs")
report_source_dir = pathlib.Path("report")
model_package_source = "model_package"
workflow_dir = "model_package/workflows"
DNS_Abaqus_dir = "model_package/DNS_Abaqus"
DNS_Ratel_dir = "model_package/DNS_Ratel"
DNS_GEOS_dir = "model_package/DNS_GEOS"
mesh_templates_dir = "model_package/meshes"
filter_source_dir = "model_package/Filter"
calibrate_source_dir = "model_package/Calibrate"
Tardigrade_MOOSE_source_dir = "model_package/Tardigrade_MOOSE"
peta_data_drive = pathlib.Path("/pl/active/psaap/share/SimulationResultsAndInputs")
peta_data_copy_abspath = project_dir / "peta_data_copy"
project_variables = {
    "project_configuration": str(project_configuration),
    "project_dir": str(project_dir),
    "project_name": project_name,
    "version": version,
    "author_list": author_list,
    "author_latex": author_latex,
    "documentation_pdf": f"{latex_project_name}-{version}.pdf",
    "report_pdf": f"{latex_project_name}-{version}-report.pdf",
    "documentation_abspath": str(project_dir / documentation_source_dir),
    "report_abspath": str(project_dir / report_source_dir),
    "model_package_abspath": str(project_dir / model_package_source),
    "workflow_abspath": str(project_dir / workflow_dir),
    "DNS_Abaqus_abspath": str(project_dir / DNS_Abaqus_dir),
    "DNS_Ratel_abspath": str(project_dir / DNS_Ratel_dir),
    "DNS_GEOS_abspath": str(project_dir / DNS_GEOS_dir),
    "mesh_templates_abspath": str(project_dir / mesh_templates_dir),
    "filter_source_abspath": str(project_dir / filter_source_dir),
    "calibrate_source_abspath": str(project_dir / calibrate_source_dir),
    "Tardigrade_MOOSE_source_abspath": str(project_dir / Tardigrade_MOOSE_source_dir),
    "micromorphic_abspath": str(micromorphic_location),
    "constraints_abspath": str(constraints_location),
    "peta_data_drive": peta_data_drive,
    "peta_data_copy_abspath": peta_data_copy_abspath,
    "regression_alias": "regression"
}
for key, value in project_variables.items():
    env[key] = value

# Make the model package and calibration tools importable for: 
# (1) SConscript files and (2) Python and Abaqus Python environments
sys.path.insert(0, str(project_dir))
sys.path.append(str(filter_location))
sys.path.append(str(micromorphic_location))
sys.path.append(str(project_dir / filter_source_dir))
sys.path.append(str(constraints_location))
sys.path.append(str(project_dir / model_package_source))
env.PrependENVPath("PYTHONPATH", str(filter_location))
env.PrependENVPath("PYTHONPATH", str(project_dir))
env.PrependENVPath("PYTHONPATH", str(micromorphic_location))
env.PrependENVPath("PYTHONPATH", str(project_dir / filter_source_dir))
env.PrependENVPath("PYTHONPATH", str(constraints_location))
env.PrependENVPath("PYTHONPATH", str(project_dir / model_package_source))

# Build path object for extension and re-use
variant_dir_base = pathlib.Path(env["variant_dir_base"])

# Add WAVES builders
env.Append(BUILDERS={
    "AbaqusJournal": waves.scons_extensions.abaqus_journal(program=env["abaqus"]),
    "AbaqusSolver": waves.scons_extensions.abaqus_solver(program=env["abaqus"]),
    "AbaqusExtract": waves.scons_extensions.abaqus_extract(program=env["abaqus"]),
    "PythonScript": waves.scons_extensions.python_script(),
    "CondaEnvironment": waves.scons_extensions.conda_environment(),
    "RatelSolver": ratel_builder_select(),
    "TardigradeSolver": tardigrade_builder_select(),
    #"ParaviewImage": paraview_image,
    "SphinxBuild": waves.scons_extensions.sphinx_build(program=env["sphinx_build"], options="-W"),
    "SphinxPDF": waves.scons_extensions.sphinx_latexpdf(program=env["sphinx_build"], options="-W")
})
env.Append(SCANNERS=waves.scons_extensions.sphinx_scanner())

# Dump the Conda environment as documentation of as-built target environment
environment_target = env.CondaEnvironment(
    target=["environment.yaml"],
    source=[])
env.AlwaysBuild(environment_target)
Default(environment_target)

# Add documentation target(s)
if not env["ignore_documentation"]:
    # Project documentation
    build_dir = variant_dir_base / documentation_source_dir
    docs_aliases = SConscript(str(documentation_source_dir / "SConscript"), variant_dir=str(build_dir), exports=["env", "project_variables"])

    # Analysis report
    report_dir = pathlib.Path("report")
    build_dir = variant_dir_base / report_dir
    SConscript(str(report_dir / "SConscript"),
               variant_dir=str(build_dir),
               exports=["env", "project_variables"],
               duplicate=True)
else:
    print(f"The 'ignore_documentation' option was set to 'True'. Skipping documentation SConscript file(s)")
    docs_aliases = []

# Add simulation targets
workflow_configurations = [
    # Abaqus quasistatic elastic cylinder
    "Abaqus_elastic_cylinder",
    "Abaqus_elastic_cylinder_multi_domain",
    "Abaqus_elastic_cylinder_clamped",
    "Abaqus_elastic_cylinder_clamped_multi_domain",
    # Abaqus dynamic implicit elastic cylinder
    "Abaqus_elastic_cylinder_dynamic_imp",
    "Abaqus_elastic_cylinder_dynamic_imp_multi_domain",
    # Ratel quasistatic elastic cylinder
    "Ratel_elastic_cylinder",
    "Ratel_elastic_cylinder_multi_domain",
    "Ratel_elastic_cylinder_clamped",
    "Ratel_elastic_cylinder_clamped_multi_domain",
    # Ratel F83 workflows
    "Ratel_F83_multi_domain",
    # Ratel I41_02 workflows
    "Ratel_I41_02_elastic_multi_domain",
    "Ratel_I41_02_elastic_single_RVEs",
    # Tardigrade solo studies
    "Tardigrade_convergence",
]
for workflow in workflow_configurations:
    build_dir = str(variant_dir_base / workflow)
    workflow_sconscript = pathlib.Path(f"{workflow_dir}/{workflow}")
    SConscript(workflow_sconscript, variant_dir=build_dir, exports="env", duplicate=False)

# Update local copies of Peta data
import model_package.peta
if env["peta_data_copy"]:
    try:
        model_package.peta.peta_copy(peta_data_drive, peta_data_copy_abspath)
    except Exception as err:
        print("Peta data copy failed!")

# Configure softare
import model_package.config_software
if env["config_software"]:
    try:
        model_package.config_software.config_software(config_file)
    except Exception as err:
        print("Software configuration failed!")

# Project title
print(r"""\
█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗
╚════╝╚════╝╚════╝╚════╝╚════╝╚════╝╚════╝╚════╝╚════╝╚════╝╚════╝╚════╝╚════╝╚════╝


    ████████╗ █████╗ ██████╗ ██████╗ ██╗ ██████╗ ██████╗  █████╗ ██████╗ ███████╗
    ╚══██╔══╝██╔══██╗██╔══██╗██╔══██╗██║██╔════╝ ██╔══██╗██╔══██╗██╔══██╗██╔════╝
       ██║   ███████║██████╔╝██║  ██║██║██║  ███╗██████╔╝███████║██║  ██║█████╗
       ██║   ██╔══██║██╔══██╗██║  ██║██║██║   ██║██╔══██╗██╔══██║██║  ██║██╔══╝
       ██║   ██║  ██║██║  ██║██████╔╝██║╚██████╔╝██║  ██║██║  ██║██████╔╝███████╗
       ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝ ╚══════╝


                         ██████████████████████████████
                    ███████      █████               ██████████
                 ████                 █                █      █████
               ███                     █                █          ████
              ██                        █                █            ███
             ██                          █                █     █       ██
    ██████████ █                          █                █     █ █     ███
  ███      █    █                          █               ██     █ █      █
  ██        █    █                         █                █       ██     ██
 ██         ██   ██                        ██               ██       █      ██
  █         █    █                         █               ██       █    █   ██
  ██       █    █                          █               █       █    █  █  ██
   ██     █    █                          █              ██       █    █  █    ██
    ██████████                        ██ ██              █    █  █    █  █     ██
             ██                     ██  ██              ██      ██    █ █      ██
              ███                ███   ██              ██    █ █      █ █      ██
                ████   ██████████    ███            █ ██     █        █        ██
                  ██████                            █ █        █              ██
                     █████████     ███              █ █      █ █     █ █      ██
                      █      ██       ██              █      █ █     █ █      █
                      █       ██       ██           █ █      █       █ █      █
                     ██        ███      █             █ █         █    █     ██
                      ██        █████████            █ █      █ █     █     ██
                       ██   ████        █            █ █            █ ███████
                       ██████           ██            █      ███   ██  █ █ █
                        █ █ █            ██           █      █ ██████
                                          ██         ███    ██  █ █ █
                                           █       ██   █████
                                            ███████     █ █ █
                                             █ █ █


        ███████╗██╗  ██╗ █████╗ ███╗   ███╗██████╗ ██╗     ███████╗███████╗
        ██╔════╝╚██╗██╔╝██╔══██╗████╗ ████║██╔══██╗██║     ██╔════╝██╔════╝
        █████╗   ╚███╔╝ ███████║██╔████╔██║██████╔╝██║     █████╗  ███████╗
        ██╔══╝   ██╔██╗ ██╔══██║██║╚██╔╝██║██╔═══╝ ██║     ██╔══╝  ╚════██║
        ███████╗██╔╝ ██╗██║  ██║██║ ╚═╝ ██║██║     ███████╗███████╗███████║
        ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝     ╚══════╝╚══════╝╚══════╝


█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗█████╗
""")

# Add default target list to help message
# Add aliases to help message so users know what build target options are available
# This must come *after* all expected Alias definitions and SConscript files.
waves.scons_extensions.project_help_message()

# Write targets to .scon_autocomplete file
with open('.scons_autocomplete', 'w') as f:
    f.write(' '.join(default_ans))