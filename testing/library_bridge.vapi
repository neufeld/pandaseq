[CCode(instance_pos = 0.1, cname = "Assemble", cheader_filename = "library_bridge.h")]
public extern delegate unowned Panda.result_seq? Assemble(Panda.identifier id, [CCode(array_length_type = "size_t")]Panda.qual[] forward, [CCode(array_length_type = "size_t")]Panda.qual[] reverse);
