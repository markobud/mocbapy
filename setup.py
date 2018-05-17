from setuptools import setup, Extension

setup(
    name="mocbapy",
    version='0.2',
    description='Multi Objective Constrait-Based Analysis in python',
    author='Marko Budinich',
    author_email='marko.budinich@ls2n.fr',
    url='https://gitlab.univ-nantes.fr/mbudinich/mocbapy',
    license='GPLv3',
    long_description='Multi Objective Constrait-Based Analysis in python',
    packages=['mocbapy'],
    platforms='linux, OSX', install_requires=['pandas', 'optlang', 'cobra', 'tqdm', 'numpy', 'scipy']
)
