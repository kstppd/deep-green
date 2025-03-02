project(
    'eulernv',
    ['cuda', 'cpp'],
    version: '3.0',
    default_options: ['warning_level=3', 'cpp_std=c++20'],
)

CUSTOM_INSTALL_PATH='./external/'
runner = run_command('git', 'clone', 'https://github.com/highfive-devs/highfive.git', 'external/highfive', check: false)
spdlog_dep = dependency('spdlog')
cuda_dep = dependency('cuda', version: '>=11', modules: ['nvToolsExt'])
hdf5_dep = dependency('hdf5', language: 'cpp')
includes=[CUSTOM_INSTALL_PATH+'./highfive/include/']
deps = [cuda_dep, hdf5_dep,spdlog_dep]

add_project_arguments(
    [
        '--expt-relaxed-constexpr',
        '-std=c++20',
        '--expt-extended-lambda',
        '-Xcompiler',
        '-fopenmp',
        '-Xcompiler',
        '-Wno-pedantic',
        '--std=c++20',
        '-Xcudafe',
        '--diag_suppress=20208',
    ],
    language: 'cuda',
)

executable(
    'eulernv',
    sources: ['eulernv.cu'],
    include_directories: includes,
    link_args: '-fopenmp',
    dependencies: deps,
    install: true,
)
