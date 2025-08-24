import sys
def write_data(N_ROD, N_SEG, L, H, K3, LOWER_RAT, F_I, F_S, BE, RAND_MAG):
    with open("data.in", "w") as f:
        print("in?")
        f.write(f"{N_ROD}\n")
        f.write(f"{N_SEG}\n")
        f.write(f"{L}\n")
        f.write(f"{H}\n") 
        f.write(f"{K3}\n")
        f.write(f"{LOWER_RAT}\n")
        f.write(f"{F_I}\n")
        f.write(f"{F_S}\n")
        f.write(f"{BE}\n")
        f.write(f"{RAND_MAG}\n")
        f.close()       

write_data(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])