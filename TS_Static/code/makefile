GOAL = dip_10000.x

INCL = common.inc

OBJ = boxreset.o\
main.o\
kset.o\
conjrod.o\
dcell.o\
dipolemove.o\
dofmoveselect.o\
dump.o\
ener.o\
erfunc.o\
gauss.o\
header.o\
invert.o\
matinv.o\
output.o\
tramp.o\
pramp.o\
prescalcbo.o\
ran1.o\
rattle.o\
rdftern.o\
rdfouttern.o\
readin.o\
realE_multipol_ewald_nosr.o\
recipE_tern_ewald_varkmax.o\
sr_energy.o\
rescale.o\
rescalen.o\
rgdrealE_ewald_nosr.o\
rgdrecipE_tern_ewald.o\
rstrun.o\
setup.o\
shake.o\
debye_scherer.o\
debye_scherer_out.o\
timeset.o\
trans_chains.o\
update_positions_fixed.o\
update_velocities_fixed.o\
velkill.o\
vellkill.o\
velsetvv.o

#FFLAGS = -O2 -static
# FFLAGS = -Wall
 FFLAGS = -O
#LIBS = /home/coulson/mw/mwilson/libs/libnaglin.a
                                                                                        
${GOAL}: ${OBJ} ; \
   gfortran ${OBJ} ${FFLAGS} ${LIBS} -o ${GOAL}
#  f77 ${OBJ} ${FFLAGS} ${LIBS} -o ${GOAL}
#         pgf90   ${OBJ} ${FFLAGS} -o ${GOAL}
                                                                                        
.f.o: ${INCL} ; \
  gfortran  ${FFLAGS} -c $<
#  f77  ${FFLAGS} -c $<
#        pgf90   ${FFLAGS} -c $<
                                                                                        
program:;\
       cat ${SRCS} > comp
