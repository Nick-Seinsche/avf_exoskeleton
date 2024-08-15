function test_A3_sparse()
    A = assemble_A3(3);
    B = assemble_A3_old(3);

    A - B