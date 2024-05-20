project(
  'stringspinner',
  'cpp', 'fortran',
  meson_version: '>=1.2',
  default_options: {
    'cpp_std':   'c++11',
    'buildtype': 'release',
    'libdir':    'lib',
  },
  version: 'Comput.Phys.Commun. 292 (2023) 108886', # no version number, so use the publication reference instead; arXiv: https://arxiv.org/abs/2305.05058
)

# dependency pythia8
pythia_link_args = []
pythia_config_prog = find_program('pythia8-config', required: false)
if not pythia_config_prog.found()
  error('Pythia8 is not found')
endif
foreach arg : run_command(pythia_config_prog, '--libs', check: true).stdout().strip().split()
  if not arg.contains('-rpath')
    pythia_link_args += arg
  endif
endforeach
pythia_version_res = run_command(pythia_config_prog, '--version', check: false)
pythia_version = pythia_version_res.stdout().strip()
if pythia_version_res.returncode() != 0 or pythia_version == ''
  pythia_version = '8' # fallback; see https://gitlab.com/Pythia8/releases/-/issues/448
endif
stringspinner_pythia_dep = declare_dependency(
  include_directories: run_command(pythia_config_prog, '--includedir', check: true).stdout().strip(),
  link_args: pythia_link_args,
  version: pythia_version,
)

# headers
install_headers(
  [
    'stringspinner/PrimordialKT.h',
    'stringspinner/StringSpinner.h',
    'stringspinner/Transversity.h',
    'stringspinner/VectorMesonDecays.h',
  ],
  subdir: meson.project_name(),
)

# libraries
def_lib = static_library(
  'StringSpinner_def',
  'stringspinner/definitions.f90',
  install: true,
  fortran_args: ['-Wno-tabs'],
)
mc3P0_lib = static_library(
  'StringSpinner_mc3P0',
  'stringspinner/mc3P0.f90',
  install: true,
  link_with: def_lib,
  fortran_args: ['-frecord-marker=8', '-fbounds-check', '-Wno-tabs', '-Wno-unused-variable', '-Wno-conversion', '-Wno-unused-dummy-argument', '-Wno-aliasing'],
)
stringspinner_libs = [def_lib, mc3P0_lib]

# example executable
if get_option('install_example')
  dis_exe = executable(
    'stringspinner-example',
    sources: ['stringspinner/dis.cc'],
    link_with: stringspinner_libs,
    dependencies: [ pythia_dep ],
    install: true,
    cpp_args: ['-Wno-sign-compare', '-Wno-unused-variable'],
  )
endif

# declare subproject
stringspinner_dep = declare_dependency(
  include_directories: include_directories('.', is_system: true), # use `is_system` to suppress warnings in consumer builds
  link_with: stringspinner_libs,
  dependencies: [stringspinner_pythia_dep],
)
