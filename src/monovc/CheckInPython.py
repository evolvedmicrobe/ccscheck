from ConsensusCore2 import *

g = MonoMolecularIntegrator("GCAGCCTA", IntegratorConfig(), SNR(7.0, 7.0, 7.0, 7.0), "P6-C4")



mutations = [Mutation(MutationType_DELETION, 1, "C"),
             Mutation(MutationType_DELETION, 1, "G"),
             Mutation(MutationType_DELETION, 1, "A")]
print str(g)
for mut in mutations:
    g.ApplyMutation(mut)
print str(g)

for mut in mutations:
    g.ApplyMutation(Mutation(MutationType_INSERTION, 1, mut.Base))
print str(g)

