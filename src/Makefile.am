AUTOMAKE_OPTIONS = subdir-objects
AM_CPPFLAGS = @INCLUDE@ @CPPFLAGS@
AM_CXXFLAGS = @CXXFLAGS@
AM_LDFLAGS  = @LDFLAGS@

noinst_HEADERS = \
	basic_bull.hpp \
	bull.hpp \
	bull_io.hpp \
	bull_kernel.hpp \
	bull_listener.hpp \
	bull_parsers.hpp \
	char_encoding.hpp \
	complex.hpp \
	encoded_chars.hpp \
	genetic_codes.hpp \
	jph.hpp \
	like_attribute_sets.hpp \
	like_attributes.hpp \
	lin_alg.hpp \
	matrices.hpp \
	model.hpp \
	mth_exception.hpp \
	node.hpp \
	parameter.hpp \
	part_mod_index.hpp \
	settings.hpp \
	ssrf_codon_sub_mod.hpp \
	string_extensions.hpp \
	tools.hpp \
	tree.hpp \
	util.hpp \
	xbull.hpp 

bin_PROGRAMS = bull

bull_SOURCES = \
	bull.cpp \
	bull_io.cpp \
	bull_kernel.cpp \
	bull_parsers.cpp \
	char_encoding.cpp \
	complex.cpp \
	genetic_codes.cpp \
	like_attribute_sets.cpp \
	like_attributes.cpp \
	lin_alg.cpp \
	main.cpp \
	matrices.cpp \
	model.cpp \
	node.cpp \
	parameter.cpp \
	settings.cpp \
	ssrf_codon_sub_mod.cpp \
	string_extensions.cpp \
	tools.cpp \
	tree.cpp \
	util.cpp

bull_LDADD =  $(LDADD) @BULL_LIBS@

