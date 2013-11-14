[CCode(cname = "panda_assembler_new_from_file", cheader_filename = "test_wrapper.h")]
public extern Assemble? panda_assembler_new_from_file(Panda.LogProxy logger);

[CCode(instance_pos = 0.1, cname = "Assemble", cheader_filename = "test_wrapper.h")]
public extern delegate unowned Panda.result_seq? Assemble(Panda.identifier id, [CCode(array_length_type = "size_t")]Panda.qual[] forward, [CCode(array_length_type = "size_t")]Panda.qual[] reverse);


