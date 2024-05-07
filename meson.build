project(
  'clas-stringspinner',
  'cpp',
  meson_version: '>=1.2',
  default_options: {
    'cpp_std':   'c++17',
    'buildtype': 'release',
    'libdir':    'lib',
  },
  version: '0.0.0',
)

# subprojects
stringspinner_proj = subproject('stringspinner')

# dependencies
fmt_dep           = dependency('fmt', method: 'pkg-config', version: '>= 9.1.0')
pythia_dep        = dependency('pythia', fallback: ['stringspinner', 'stringspinner_pythia_dep'], version: '>=8')
stringspinner_dep = dependency('stringspinner', fallback: ['stringspinner', 'stringspinner_dep'])
