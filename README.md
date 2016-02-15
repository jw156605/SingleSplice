# SingleSplice
Algorithm for detecting alternative splicing in a population of single cells

System Requirements:
- R, Perl, gcc, and git
- Boost C++ library

Installation instructions:
- Check that your system meets the requirements
- Navigate to the desired installation directory
- Clone the SingleSplice repository:
  git clone https://github.com/jw156605/SingleSplice
- To install the Boost library, you can run the following commands:
  cd SingleSplice/
  wget https://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz
  tar -xvf boost_1_60_0.tar.gz
- If you install in a directory other than SingleSplice/boost_1_60_0/, you need to modify the BOOST
  variable in the Makefile to point to the installation directory.
- Run the Makefile by simply typing: make