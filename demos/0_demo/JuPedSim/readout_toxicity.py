import matplotlib.pyplot as plt
import numpy as np
import xml.etree.ElementTree as ET
import math

def label(quantity):
	if quantity == 'E':
		return 'Extinction Coefficient (1/m)'
	if quantity == 'FED_In':
		return 'FED Incapacitation'

	#todo:
	#"CO2 concentration (ppm)","CO concentration (ppm)",
	#"HCN concentration (ppm)", "HCl concentration (ppm)", "FED Heat Dose",
	#"FEC Smoke", "FIC Impairment", "FIC Incapacitation"

########## change here ##########
#number of agents
N = 50 #todo: parse from traj.xml
# possible quantities are: 'E','c_CO2','c_CO','c_HCl','c_HCN','FED_In','FEC_Smoke','FED_Heat','FIC_Im','FED_In'
quantity_select = 'E'
#################################

### traj XML readin
tree = ET.parse('toxicity_output_jpsfire.xml')
root = tree.getroot()

### readout
data = []
for frame in root.iter('frame'):
    frame_id = int(frame.get('ID'))
    for agent in frame.iter('agent'):
    	agent_id = int(agent.get('ID'))
    	t = float(agent.get('t'))
    	quantity = float(agent.get('%s' % quantity_select))
    	data += [agent_id, t, quantity]
data = np.array(data).reshape((-1, 3))

for i in range(1,N+1):
	
	ids_all = data[:,0]
	time_all = data[:,1]
	quantity_all = data[:,2]

	time_max = np.max(time_all)
	quantity_max = np.max(quantity_all)

	if i == 1:
		plt.plot(time_all[ids_all==1], quantity_all[ids_all==1], 'b', lw=2, label = 'agent time series')
	plt.plot(time_all[ids_all==i], quantity_all[ids_all==i], 'b', lw=2, zorder=1)
#plt.plot([0,time_max], [1, 1], "r--", lw=2.5, label='%s threshold' %quantity_select) #threshold
plt.plot([0,time_max], [quantity_max, quantity_max], "k--", lw=2.5, label=r'%s max.' %quantity_select)
plt.xlabel("time [s]")
plt.ylabel("%s" % label(quantity_select))
plt.xlim(50) # mean premovement time
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylim()
plt.legend(loc=0)
plt.savefig('toxicity_output_plot_%s.png' % quantity_select, dpi=400)
plt.show()