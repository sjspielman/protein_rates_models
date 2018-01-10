# Save to raxml model file
with open("jc69.raxmat", "w") as raxf:
    for i in range(20):
        for j in range(20):
            if (i==j):
                rate ="0.0"
            else:
                rate = "1.0"
            raxf.write(rate + "\n")
    for i in range(20):
        raxf.write("0.05\n")
