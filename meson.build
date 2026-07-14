project(
  'clas-stringspinner',
  'cpp', 'fortran',
  version: '2.0.0',
  meson_version: '>=1.3',
  default_options: {
    'cpp_std':   'c++20',
    'buildtype': 'release',
    'libdir':    'lib',
  },
)

##################################################################################
# dependencies
##################################################################################

fmt_dep  = dependency('fmt',   method: 'pkg-config', version: '>= 9.1.0')
hipo_dep = dependency('hipo4', method: 'pkg-config', version: '>= 4.2.0', required: false)

##################################################################################
# preprocessor macros
##################################################################################

add_project_arguments('-D CLAS_STRINGSPINNER_VERSION="' + meson.project_version() + '"', language: ['cpp'])
if hipo_dep.found()
  add_project_arguments('-D CLAS_STRINGSPINNER_HIPO', language: ['cpp'])
endif

##################################################################################
# resolve pythia8 dependency
##################################################################################

# get linker args
pythia_link_args = []
pythia_config_prog = find_program('pythia8-config', required: true)
foreach arg : run_command(pythia_config_prog, '--libs', check: true).stdout().strip().split()
  if not arg.contains('-rpath')
    pythia_link_args += arg
  endif
endforeach

# get version number
pythia_version_res = run_command(pythia_config_prog, '--version', check: false)
pythia_version = pythia_version_res.stdout().strip()
if pythia_version_res.returncode() != 0 or pythia_version == ''
  warning('pythia version auto-detection failed, since it is likely older than 8.312; setting version to 8.3')
  pythia_version = '8.3' # fallback for versions < 8.312; see https://gitlab.com/Pythia8/releases/-/issues/448
endif

# declare dependency
pythia_dep = declare_dependency(
  include_directories: include_directories(
    run_command(pythia_config_prog, '--includedir', check: true).stdout().strip(),
    is_system: true, # use `is_system` to suppress warnings in consumer builds
  ),
  link_args: pythia_link_args,
  version: pythia_version,
)

# verify version number
pythia_version_min = '8.313'
if pythia_dep.version().version_compare(f'<@pythia_version_min@')
  error(f'pythia version @pythia_version@ is too old (minimum=@pythia_version_min@)')
endif

##################################################################################
# build stringspinner submodule
##################################################################################

# # apply patches
# # - `stringspinner` disallow forks, so instead we apply patches, but that leaves the submodule repo dirty
# # - currently there are no patches, but the code is commented out here in case we need to patch in the future
# fs = import('fs')
# if not fs.is_file('stringspinner' / 'PATCHED')
#   message('APPLYING PATCHES TO SUBMODULE "stringspinner"')
#   run_command(
#     'patch', '-p1',
#     '-d', meson.project_source_root() / 'stringspinner',
#     '-i', meson.project_source_root() / 'patches' / 'canVetoFragmentation.patch',
#     check: true,
#   )
#   run_command(
#     'touch',
#     meson.project_source_root() / 'stringspinner' / 'PATCHED',
#     check: true,
#   )
# endif

# headers
install_headers(
  [
    'stringspinner' / 'CorrelationCoefficients.h',
    'stringspinner' / 'PrimordialKT.h',
    'stringspinner' / 'StringSpinner.h',
    'stringspinner' / 'Transversity.h',
    'stringspinner' / 'VectorMesonDecays.h',
  ],
  subdir: 'stringspinner'
)

# libraries
stringspinner_def_lib = static_library(
  'StringSpinner_def',
  'stringspinner' / 'definitions.f90',
  install: true,
  fortran_args: [
    '-Wno-tabs',
  ],
)
stringspinner_mc3P0_lib = static_library(
  'StringSpinner_mc3P0',
  'stringspinner' / 'mc3P0.f90',
  install: true,
  link_with: stringspinner_def_lib,
  fortran_args: [
    '-frecord-marker=8',
    '-fbounds-check',
    '-Wno-tabs',
    '-Wno-unused-variable',
    '-Wno-conversion',
    '-Wno-unused-dummy-argument',
    '-Wno-aliasing',
  ],
)

# declare dependency
stringspinner_dep = declare_dependency(
  include_directories: include_directories(
    'stringspinner',
    is_system: true, # use `is_system` to suppress warnings in consumer builds
  ),
  link_with: [stringspinner_def_lib, stringspinner_mc3P0_lib],
  version: '2.0', # see https://arxiv.org/pdf/2606.27162
  dependencies: [pythia_dep],
)

# optionally install examples
if get_option('install_stringspinner_example')
  foreach exe_name : ['dis', 'eeToqq', 'ee']
    executable(
      f'stringspinner-example-@exe_name@',
      sources: ['stringspinner' / f'@exe_name@.cc'],
      dependencies: [ stringspinner_dep, pythia_dep ],
      install: true,
      cpp_args: [
        '-Wno-sign-compare',
        '-Wno-unused-variable',
        '-Wno-unused-but-set-variable',
        '-Wno-maybe-uninitialized',
      ],
    )
  endforeach
endif

##################################################################################
# main executable
##################################################################################

# build main executable
main_exe = executable(
  meson.project_name(),
  'src' / meson.project_name() + '.cpp',
  'src' / 'config' / 'src' / 'common.cc',
  dependencies: [ fmt_dep, stringspinner_dep, pythia_dep, hipo_dep ],
  install: true,
)

# test main executable
test(
  meson.project_name(),
  main_exe,
  args: [
    '--num-events', '10000',
    '--pol-type', 'LU',
    '--beam-spin', 'p',
    '--cut-inclusive', '11,211,-211',
    '--seed', '828203788',
  ]
)

# install scripts
install_data(
  'src' / 'clas-stringspinner-slurm.rb',
  rename: 'clas-stringspinner-slurm',
  install_dir: get_option('bindir'),
  install_mode: 'rwxr-xr-x',
)

##################################################################################
# summary
##################################################################################

summary(
  {
    'pythia':        pythia_dep.version(),
    'stringspinner': stringspinner_dep.version(),
  },
  section: 'Dependency Versions:'
)
