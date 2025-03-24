project(
    'eulernv',
    ['cuda', 'cpp'],
    version: '3.0',
    default_options: ['warning_level=3', 'cpp_std=c++20'],
)

CUSTOM_INSTALL_PATH='./external/'
runner = run_command('git', 'clone', 'https://github.com/highfive-devs/highfive.git', 'external/highfive', check: false)
runner = run_command('git', 'clone', 'https://github.com/kstppd/asterix.git', 'external/asterix', check: false)
runner = run_command('git', 'clone', 'https://github.com/kstppd/libnpy.git', 'external/libnpy', check: false)
spdlog_dep = dependency('spdlog')
cuda_dep = dependency('cuda', version: '>=11', modules: ['cublas', 'nvToolsExt'])
hdf5_dep = dependency('hdf5', language: 'cpp')
blas_dep = dependency('openblas')
includes=[CUSTOM_INSTALL_PATH+'./highfive/include/',CUSTOM_INSTALL_PATH+'./asterix/include/',CUSTOM_INSTALL_PATH+'./libnpy/include/']
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

executable(
    'neural_sdf',
    sources: ['neural_sdf.cu'],
    include_directories: includes,
    link_args: '-fopenmp',
    dependencies: [cuda_dep, spdlog_dep,blas_dep],
    install: true,
)
