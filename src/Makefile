all: occu_tmb_re.so

occu_tmb_re.so: occu_tmb_re.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('occu_tmb_re.cpp','-O0 -g -Wno-ignored-attributes', safebounds=FALSE)"

clean:
	rm -rf *o
