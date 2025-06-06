project(
  'clas-stringspinner',
  'cpp',
  meson_version: '>=1.3',
  default_options: {
    'cpp_std':    'c++20',
    'buildtype':  'release',
    'libdir':     'lib',
    'wrap_mode':  'forcefallback',
  },
  version: '0.1.5',
)

# subprojects
stringspinner_proj = subproject('stringspinner')

# dependencies
fmt_dep           = dependency('fmt', method: 'pkg-config', version: '>= 9.1.0')
stringspinner_dep = dependency('stringspinner', fallback: ['stringspinner', 'stringspinner_dep'])
pythia_dep        = dependency('pythia', fallback: ['stringspinner', 'stringspinner_pythia_dep'], version: '8.312')

# main executable
add_project_arguments('-DCLAS_STRINGSPINNER_VERSION="' + meson.project_version() + '"', language: ['cpp'])
main_exe = executable(
  meson.project_name(),
  'src' / meson.project_name() + '.cpp',
  'src' / 'config' / 'src' / 'common.cc',
  dependencies: [ fmt_dep, stringspinner_dep, pythia_dep ],
  install: true,
)

# test executable
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

###################################################################################################################
# ————————————————
# 1) C++ compiler object
# ————————————————
cxx = meson.get_compiler('cpp')

# ————————————————
# 2) HIPO4 (via pkg-config)
# ————————————————
hipo4_dep = dependency('hipo4', version: '>=4.2.0', required: true)

# ————————————————
# 3) CLAS12ROOT (hardcoded paths)
# ————————————————
clas12root_inc = include_directories(
  '/cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.5/4.2.0/Clas12Root',
  '/cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.5/4.2.0/Clas12Banks',
)

clas12root_libdir = '/cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.5/4.2.0/lib'
clas12root_lib    = cxx.find_library('Clas12Root',  dirs: [clas12root_libdir], required: true)
clas12banks_lib   = cxx.find_library('Clas12Banks', dirs: [clas12root_libdir], required: true)

clas12_dep = declare_dependency(
  include_directories: clas12root_inc,
  dependencies: [ clas12root_lib, clas12banks_lib ]
)

# ————————————————
# 4) ROOT flags (via root-config)
# ————————————————
root_cflags = run_command('root-config', '--cflags', check: true).stdout().strip().split()
root_libs   = run_command('root-config', '--libs',   check: true).stdout().strip().split()

# ————————————————
# 5) Local project includes
# ————————————————
src_inc = include_directories('src')

# ————————————————
# 6) Build the hipo2tree executable
# ————————————————
hipo2tree_exe = executable(
  'hipo2tree',
  [ 'src/hipo2tree.cpp', 'src/Kinematics.cpp' ],
  include_directories: [ src_inc, clas12root_inc ],
  dependencies:        [ hipo4_dep, clas12_dep ],
  cpp_args:            root_cflags,
  link_args:           root_libs,
  install:             true,
)

# ————————————————
# 7) Build the hipo2tree executable
# ————————————————
lund2tree_exe = executable(
  'lund2tree',
  [ 'src/lund2tree.cpp', 'src/Kinematics.cpp' ],
  include_directories: [ src_inc, clas12root_inc ],
  dependencies:        [ hipo4_dep, clas12_dep ],
  cpp_args:            root_cflags,
  link_args:           root_libs,
  install:             true,
)
###################################################################################################################
