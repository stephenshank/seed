from setuptools import setup


setup(
    name='seed',
    python_requires='>3.6.6',
    version='0.0.1',
    url='https://github.com/stephenshank/seed',
    download_url="https://github.com/stephenshank/seed/archive/v0.0.1.tar.gz",
    description='Somatic epigenetic evolution detector',
    author='Stephen D. Shank',
    author_email='sshank314@gmail.com',
    maintainer='Stephen D. Shank',
    maintainer_email='sshank314@gmail.com',
    install_requires=[
        'pysam>=0.21.0',
    ],
    packages=['seed'],
    entry_points={
        'console_scripts': [
            'seed = seed.epipolymorphism:cli'
        ]
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
    ],
    include_package_data=True
)
