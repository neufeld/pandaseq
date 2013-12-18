PANDASEQ_PREFIX=/usr/local

pandaseq-diff: pandaseq-diff.vala library_bridge.h library_bridge.vapi access_new_library.vapi access_new_library.h renamed_access_new_library.o renamed_pandaseq.a
	PKG_CONFIG_PATH=${PANDASEQ_PREFIX}/lib/pkgconfig valac --vapidir ${PANDASEQ_PREFIX}/share/vala/vapi pandaseq-diff.vala library_bridge.vapi access_new_library.vapi --pkg pandaseq-2 --pkg pandaseq-2-url -X renamed_access_new_library.o -X renamed_pandaseq.a -g --save-temps

access_new_library.c access_new_library.h access_new_library.vapi: access_new_library.vala library_bridge.h library_bridge.vapi
	PKG_CONFIG_PATH=${PANDASEQ_PREFIX}/lib/pkgconfig valac --vapidir ${PANDASEQ_PREFIX}/share/vala/vapi access_new_library.vala library_bridge.vapi --pkg pandaseq-2 -g -C -H access_new_library.h --vapi access_new_library.vapi

access_new_library.o: access_new_library.c access_new_library.h library_bridge.h
	gcc -c -o $@ $< $(shell PKG_CONFIG_PATH=${PANDASEQ_PREFIX}/lib/pkgconfig pkg-config --cflags glib-2.0 pandaseq-2)

renamed_pandaseq.a: ../.libs/libpandaseq.a
	nm $< --defined-only -B -g | awk 'NF == 3 { print $$3, "xx" $$3 }' | objcopy --redefine-syms=/dev/stdin $< $@

renamed_access_new_library.o: access_new_library.o
	nm $< -B -g | awk 'NF == 2 && $$2 ~ /^panda_/ { print $$2, "xx" $$2 }' | objcopy --redefine-syms=/dev/stdin $< $@

clean:
	rm -f pandaseq-diff *.o *.a access_new_library.c access_new_library.h access_new_library.vapi

.PHONY: clean
