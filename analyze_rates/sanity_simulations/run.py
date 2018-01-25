"""
    We simulate a quick sanity with LG and JC69 to compare variance in rate inference
"""

from pyvolve import *
from phyphy import *
treestring = "(((t112:0.092,((t7:0.044,t147:0.055):0.019,((((t145:0.090,t103:0.055):0.025,t150:0.046):0.017,(t130:0.053,t65:0.082):0.065):0.089,((t49:0.060,((t80:0.086,t118:0.061):0.091,t105:0.006):0.057):0.072,((t128:0.095,(((t34:0.054,t4:0.099):0.011,t104:0.030):0.007,t146:0.034):0.035):0.015,(t27:0.029,(((t109:0.067,t76:0.068):0.007,t119:0.034):0.012,t28:0.083):0.037):0.049):0.001):0.053):0.065):0.098):0.071,(t45:0.048,t115:0.074):0.057):0.011,(((((t107:0.042,((t8:0.063,t42:0.077):0.057,t135:0.081):0.007):0.083,(((((t68:0.054,t55:0.022):0.045,(((t120:0.002,t114:0.053):0.000,t47:0.090):0.082,(t124:0.038,(t126:0.098,t5:0.014):0.035):0.089):0.097):0.030,(t44:0.046,(t111:0.069,t64:0.087):0.072):0.083):0.038,(t73:0.012,(t101:0.034,(t40:0.006,t116:0.007):0.083):0.054):0.063):0.029,((t88:0.040,(t2:0.024,(t61:0.099,t9:0.021):0.033):0.018):0.027,(((((t46:0.031,(((t48:0.010,t16:0.095):0.028,t137:0.014):0.004,t43:0.085):0.095):0.039,t93:0.009):0.071,t70:0.058):0.072,((((((t148:0.058,t79:0.045):0.057,t83:0.091):0.008,t75:0.044):0.073,(t140:0.065,t113:0.048):0.088):0.038,(((t18:0.027,t84:0.050):0.048,(t6:0.084,t26:0.091):0.038):0.071,t138:0.070):0.065):0.079,(t23:0.085,t50:0.005):0.022):0.012):0.026,(((t17:0.077,t24:0.013):0.003,t89:0.020):0.027,t97:0.074):0.034):0.047):0.024):0.016):0.032,((((((((((t39:0.079,t32:0.099):0.074,t29:0.084):0.026,(t60:0.002,t127:0.094):0.034):0.081,(t91:0.079,((t33:0.009,t98:0.096):0.065,t66:0.002):0.086):0.004):0.033,t53:0.069):0.021,((((t52:0.081,t31:0.040):0.078,t106:0.097):0.025,(t3:0.080,t121:0.025):0.054):0.016,(((t62:0.083,(t81:0.039,(t72:0.035,t71:0.050):0.054):0.083):0.003,t90:0.032):0.003,t108:0.008):0.028):0.096):0.032,(t11:0.018,((t85:0.049,(t77:0.052,t14:0.041):0.083):0.070,(t131:0.019,t35:0.048):0.062):0.066):0.065):0.028,(((t86:0.083,((t134:0.028,t37:0.066):0.087,(t21:0.020,t125:0.093):0.002):0.077):0.075,(((t133:0.056,t117:0.006):0.093,t102:0.046):0.065,t136:0.083):0.050):0.017,((t41:0.045,t122:0.087):0.034,((t144:0.030,t143:0.012):0.020,t38:0.025):0.048):0.043):0.078):0.043,((((t1:0.009,t99:0.050):0.018,t129:0.088):0.080,((t95:0.028,t123:0.027):0.098,(((t58:0.098,t12:0.006):0.078,t69:0.034):0.023,(t92:0.048,t13:0.003):0.071):0.075):0.090):0.046,((t59:0.096,(t56:0.019,((t51:0.095,t63:0.017):0.073,t67:0.008):0.049):0.064):0.045,(t132:0.063,t20:0.043):0.072):0.045):0.025):0.083,(((t139:0.059,(t19:0.090,((t96:0.037,t54:0.071):0.001,t142:0.069):0.047):0.082):0.061,t141:0.035):0.047,(((t36:0.074,t15:0.090):0.071,t25:0.066):0.036,(t57:0.079,t94:0.040):0.074):0.053):0.071):0.084):0.086,t22:0.055):0.004,((((t149:0.016,(t10:0.077,t78:0.047):0.023):0.055,t74:0.050):0.019,t87:0.028):0.009,((t100:0.055,(t82:0.093,t110:0.070):0.090):0.072,t30:0.066):0.074):0.064):0.098);"


tree = read_tree(tree = treestring)

for simmodel in ["JC69", "LG"]:
    print simmodel
    if simmodel == "LG":
        infmodel = "JC69"
    else:
        infmodel = "LG"

    for i in range(10):
        print "    ", i
#         m = Model(model)
#         p = Partition(models = m, size = 500)
#         e = Evolver(tree = tree, partitions = p)
#         e(ratefile = None, infofile = None, seqfile = None)
#         seqs = e.get_sequences()    
#     
        outfile = simmodel + str(i) + ".fna"
#         with open(outfile, "w") as f:
#             for id in seqs:
#                 f.write(">" + id + "\n" + seqs[id] + "\n")
#             f.write("\n" + treestring)

        h = HyPhy(suppress_log = True, quiet=True, executable = "HYPHYMPI", mpi_options = " -np 3" )
        l = LEISR(hyphy = h, data = outfile, type = "protein", model = infmodel)
        l.run_analysis()
        x = Extractor(l)
        x.extract_csv("sim_" + simmodel + "_leisr_" + infmodel + "_" + str(i) + ".csv")
   
    