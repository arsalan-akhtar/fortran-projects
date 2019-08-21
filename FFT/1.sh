#ls
rm *.o *.mod
#cd CML/
#make clean
make
#cd ..
#-O2 -fPIC -ftree-vectorize
#gfortran -c -O2 -fPIC -ftree-vectorize mpi_interfaces.F
#gfortran -c -O2 -fPIC -ftree-vectorize mpi_siesta.F90
#-----------------------------------------------------------------------
gfortran -c -O2 -fPIC -ftree-vectorize precision.F 
#-----------------------------------------------------------------------
gfortran -c -O2 -fPIC -ftree-vectorize parallel.F 
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize    m_wxml_error.f90
#gfortran -c -O2 -fPIC -ftree-vectorize    m_wxml_escape.f90
#gfortran -c -O2 -fPIC -ftree-vectorize    m_wxml_buffer.f90
#gfortran -c -O2 -fPIC -ftree-vectorize    m_wxml_array_str.f90
#gfortran -c -O2 -fPIC -ftree-vectorize    m_wxml_dictionary.f90
#gfortran -c -O2 -fPIC -ftree-vectorize    m_wxml_elstack.f90
#gfortran -c -O2 -fPIC -ftree-vectorize  -DFC_HAVE_ABORT  m_wxml_text.F90 
#gfortran -c -O2 -fPIC -ftree-vectorize    m_wxml_core.f90
#gfortran -c -O2 -fPIC -ftree-vectorize    m_wxml_overloads.f90
#gfortran -c -O2 -fPIC -ftree-vectorize    flib_wxml.f90
#gfortran -c -O2 -fPIC -ftree-vectorize    flib_wstml.f90
#gfortran -c -O2 -fPIC -ftree-vectorize    m_wcml_coma.f90
#gfortran -c -O2 -fPIC -ftree-vectorize    flib_wcml.f90


gfortran -c -O2 -fPIC -ftree-vectorize siesta_cml.f90 
gfortran -c -O2 -fPIC -ftree-vectorize sys.F 
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize m_io.f 
#gfortran -c -O2 -fPIC -ftree-vectorize alloc.F90 
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize cellsubs.f 
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize spatial.F

#gfortran -c -O2 -fPIC -ftree-vectorize printmatrix.F 
#gfortran -c -O2 -fPIC -ftree-vectorize schecomm.F
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize class_OrbitalDistribution.F90
#gfortran -c -O2 -fPIC -ftree-vectorize class_Sparsity.F90
#gfortran -c -O2 -fPIC -ftree-vectorize class_Data1D.F90
#gfortran -c -O2 -fPIC -ftree-vectorize class_SpData1D.F90 
#gfortran -c -O2 -fPIC -ftree-vectorize class_Data2D.F90 
#gfortran -c -O2 -fPIC -ftree-vectorize class_SpData2D.F90 
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize class_Geometry.F90 
#gfortran -c -O2 -fPIC -ftree-vectorize class_Pair_Geometry_SpData2D.F90 
#gfortran -c -O2 -fPIC -ftree-vectorize class_Fstack_Pair_Geometry_SpData2D.F90 
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize sparse_matrices.F 
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize domain_decom.F
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize parallelsubs.F
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize mesh.F 
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize fft1d.F 
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize m_walltime.f90 
#gfortran -c -O2 -fPIC -ftree-vectorize moreParallelSubs.F90 
#gfortran -c -O2 -fPIC -ftree-vectorize m_timer.F90 
#-----------------------------------------------------------------------
#gfortran -c -O2 -fPIC -ftree-vectorize fft.F  
#-----------------------------------------------------------------------
#FFT 
#-----------------------------------------------------------------------
gfortran -c -O2 -fPIC -ftree-vectorize defect.f90
#gfortran -c defect.f90
#gfortran parallel.o siesta_cml.o defect.o precision.o sys.o alloc.o fft.o cellsubs.o -o defect.exe -O2 -fPIC -ftree-vectorize


#-----------------------------------------------------------------------
#Dependence 
#-----------------------------------------------------------------------
<< --MULTILINE-COMMENT--
'//*
siesta_cmlsubs.o: files.o m_uuid.o parallel.o siesta_cml.o timestamp.o
sys.o: parallel.o siesta_cml.o
m_io.o: sys.o
alloc.o: m_io.o parallel.o precision.o sys.o
cellsubs.o: precision.o
spatial.o: precision.o
printmat.o: printmatrix.o # obsolete
schecomm.o: alloc.o
class_OrbitalDistribution.o: basic_type.inc
class_Sparsity.o: basic_type.inc
class_Sparsity.o: alloc.o
class_Data1D.o: class_Data1D.T90 basic_type.inc
class_Data1D.o: alloc.o
class_SpData1D.o: class_SpData1D.T90 basic_type.inc
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_OrbitalDistribution.o
class_SpData1D.o: class_Sparsity.o
class_Data2D.o: class_Data2D.T90 basic_type.inc
class_Data2D.o: alloc.o
class_SpData2D.o: class_SpData2D.T90 basic_type.inc
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o class_Data2D.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_OrbitalDistribution.o
class_SpData2D.o: class_Sparsity.o
class_Geometry.o: basic_type.inc
class_Geometry.o: alloc.o
class_Pair_Geometry_SpData2D.o: Pair.T90 basic_type.inc
class_Pair_Geometry_SpData2D.o: class_Geometry.o class_SpData2D.o
class_Fstack_Pair_Geometry_SpData2D.o: Fstack.T90 basic_type.inc
class_Fstack_Pair_Geometry_SpData2D.o: class_Pair_Geometry_SpData2D.o
sparse_matrices.o: alloc.o class_Fstack_Pair_Geometry_SpData2D.o
sparse_matrices.o: class_OrbitalDistribution.o class_SpData1D.o
sparse_matrices.o: class_SpData2D.o class_Sparsity.o precision.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o sys.o
parallelsubs.o: domain_decom.o parallel.o spatial.o sys.o
mesh.o: precision.o
fft1d.o: parallel.o precision.o sys.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o sys.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o sys.o
fft.o: alloc.o fft1d.o m_timer.o mesh.o parallel.o parallelsubs.o precision.o
fft.o: sys.o
*//'
--MULTILINE-COMMENT--
