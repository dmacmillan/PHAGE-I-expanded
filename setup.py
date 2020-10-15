from setuptools import setup, find_packages

setup(name='phagei',
      author='Daniel MacMillan',
      author_email='drm5@sfu.ca',
      version='v1.0.0',
      packages=find_packages(),
      install_requires=['pytest', 'pyyaml'],
      python_requires='>=3.6',
      package_data={'phagei': ['epitopes_v1.0.3.txt']},
      entry_points={'console_scripts': ['phagei = phagei.main:main']})
