import setuptools
from os import path
import spci

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="spci",
    version=spci.__version__,
    author="Pavel Polishchuk",
    author_email="pavel_polishchuk@ukr.net",
    description="SPCI: structural and physicochemical interpretation of QSAR models",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DrrDom/spci",
    packages=['spci'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    python_requires='>=3.6',
    install_requires=['sirms>=1.2.0'],
    extras_require={
        'rdkit': ['rdkit>=2017.09'],
    },
    entry_points={'console_scripts':
                  ['spci = spci.spci:main',
                   'spci_add_relative_frag_size = spci.utils.add_relative_frag_size:main',
                   'spci_calc_atomic_properties_chemaxon = spci.calc_atomic_properties_chemaxon:main',
                   'spci_calc_frag_contrib = spci.calc_frag_contrib:main',
                   'spci_descriptors = spci.descriptors:entry_point',
                   'spci_find_frags_auto_rdkit = spci.find_frags_auto_rdkit:main',
                   'spci_find_frags_rdkit = spci.find_frags_rdkit:main',
                   'spci_find_murcko_rdkit = spci.find_murcko_rdkit:main',
                   'spci_model = spci.model:main',
                   'spci_plot_contributions = spci.plot_contributions:main',
                   'spci_plot_property_distribution = spci.plot_property_distribution:main',
                   'spci_predict = spci.predict:main']},
    include_package_data=True
)
