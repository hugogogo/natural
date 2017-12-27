## Test environments
* local OS X 10.12.6; R 3.2.4
* CentOS 6.2; R 3.3.0
* win-builder (x86_64-w64-mingw32 (64-bit); R Under development (unstable) (2017-09-12 r73242))

## R CMD check results (OS X 10.12.6)
0 errors | 0 warnings | 1 note
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Guo Yu <gy63@cornell.edu>'

## R CMD check results (CentOS 6.2)
0 errors | 1 warnings | 1 note
* WARNING
‘qpdf’ is needed for checks on size reduction of PDFs
* checking compiled code ... NOTE
File ‘natural/libs/natural.so’:
  Found no calls to: ‘R_registerRoutines’, ‘R_useDynamicSymbols’

## R CMD check results (x86_64-w64-mingw32 (64-bit))
0 errors | 0 warnings | 2 note
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Guo Yu <gy63@cornell.edu>'

* checking compiled code ... NOTE
File 'natural/libs/i386/natural.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
File 'natural/libs/x64/natural.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
