from distutils.core import setup 
import sys

packagedata={'Pfeature': ['/*']}

setup(name='Pfeature',

	version='1.0',
	
	description="Computation of the features of the prorein and peptide sequences" ,

        long_description="Pfeature is a comprehensive package which will allow users to compute most of the protein features that have been discovered over the past decades. Using different functional modules of this package, users will be able to evaluate five major categories of protein features: i) Composition-based features, ii) Binary profile of sequences, iii) evolutionary information based features, iv) structural descriptors,and v) pattern based descriptors, for a group of protein/peptide sequences. Additionally, users will also be able to generate these features for sub-parts of protein/peptide sequences. This will be helpful to annotate structure, function and therapeutic properties of proteins.",
	
	
	author='Raghava Group',
	
	author_email='raghava@iiitd.ac.in',
	
	license='GNU General Public License v3.0',
	
	packages=['Pfeature'],
	
	zip_safe=False,
	
	package_data=packagedata,

#data_files = datafiles,

	package_dir={'Pfeature':'Pfeature'},

	scripts = [],

	py_modules = []
      
      )
      
