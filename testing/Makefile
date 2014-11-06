PANDASEQ_PREFIX=/usr/local

reg-test: reg-test.vala control_setup.o experiment_setup.o renamed_pandaseq.a
	PKG_CONFIG_PATH=${PANDASEQ_PREFIX}/lib/pkgconfig valac --vapidir ${PANDASEQ_PREFIX}/share/vala/vapi reg-test.vala --pkg pandaseq-2 --pkg pandaseq-2-url -X control_setup.o -X experiment_setup.o -X renamed_pandaseq.a -g --save-temps -X -lpthread -X -lbz2 -X -lz -X -lm -X -lltdl

setup.c: setup.vala
	PKG_CONFIG_PATH=${PANDASEQ_PREFIX}/lib/pkgconfig valac -g -C setup.vala --pkg pandaseq-2

control_setup.o: setup.c
	gcc -Dcreate_assembler=create_assembler_control -c -o $@ $< $(shell PKG_CONFIG_PATH=${PANDASEQ_PREFIX}/lib/pkgconfig pkg-config --cflags glib-2.0 pandaseq-2)

experiment_raw_setup.o: setup.c
	gcc -Dcreate_assembler=create_assembler_experiment -c -o $@ $< $(shell PKG_CONFIG_PATH=${PANDASEQ_PREFIX}/lib/pkgconfig pkg-config --cflags glib-2.0) -I..

renamed_pandaseq.a: ../.libs/libpandaseq.a
	nm $< --defined-only -B -g | awk 'NF == 3 { print $$3, "xx" $$3 }' | objcopy --redefine-syms=/dev/stdin $< $@

experiment_setup.o: experiment_raw_setup.o
	nm $< -B -g | awk 'NF == 2 && $$2 ~ /^panda_/ { print $$2, "xx" $$2 }' | objcopy --redefine-syms=/dev/stdin $< $@

clean:
	rm -f reg-test *.o *.a setup.c reg-test.c

.PHONY: clean
