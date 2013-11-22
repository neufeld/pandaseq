PANDASEQ_PREFIX=/usr/local

pandaseq-diff: pandaseq-diff.vala test_wrapper.h test_wrapper.vapi test_wrapper.c test_pandaseq.a
	PKG_CONFIG_PATH=${PANDASEQ_PREFIX}/lib/pkgconfig valac --vapidir ${PANDASEQ_PREFIX}/share/vala/vapi pandaseq-diff.vala test_wrapper.vapi test_wrapper.c --pkg pandaseq-2 --pkg pandaseq-2-url -X test_pandaseq.a -g

test_pandaseq.a: ../.libs/libpandaseq.a
	nm $< --defined-only -B -g | awk 'NF == 3 { print $$3, "xx" $$3 }' | objcopy --redefine-syms=/dev/stdin $< $@
