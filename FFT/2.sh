#gfortran defect.o precision.o parallel.o siesta_cml.o  sys.o libwxml.a -o defect  #-O2 -fPIC -ftree-vectorize
#gfortran defect.o cellsubs.o -o defect
#gfortran m_wxml_error.o m_wxml_escape.o m_wxml_buffer.o \
#m_wxml_array_str.o m_wxml_dictionary.o m_wxml_elstack.o m_wxml_text.o \
#m_wxml_core.o m_wxml_overloads.o flib_wxml.o flib_wstml.o m_wcml_coma.o \
#flib_wcml.o siesta_cml.o -o sys

#gfortran -o defect \
#libwxml.a m_wxml_buffer.o m_wxml_array_str.o m_wxml_dictionary.o \
#m_wxml_elstack.o m_wxml_text.o m_wxml_escape.o  m_wxml_core.o \
#m_wxml_overloads.o flib_wxml.o flib_wstml.o flib_wcml.o m_wxml_error.o \
#m_wcml_coma.o sys.o defect.o 



gfortran -o defect \
 alloc.o  precision.o cellsubs.o parallel.o m_io.o siesta_cml.o \
 sys.o pxf.o m_io_s.o sockets.o fsockets.o defect.o libfdf.a libwxml.a libxmlparser.a
