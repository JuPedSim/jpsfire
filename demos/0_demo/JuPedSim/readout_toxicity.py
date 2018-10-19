import matplotlib.pyplot as plt
import numpy as np
import xml.etree.ElementTree as ET
import math

### traj XML readin
tree = ET.parse('toxicity_output_jpsfire.xml')
root = tree.getroot()

### empty arrays
t = E = co2 = co = hcn = hcl = fed_in = fec_smoke = fed_heat = fic_im = fic_in = np.array([])

### readout
for frame in root.iter('frame'):
    frame_id = int(frame.get('ID'))
    for agent in frame.iter('agent'):
    	agent_id = int(agent.get('ID')) # todo: is it necessary to plot agent specifically
    	t = np.append(t, float(agent.get('t')))
    	E = np.append(E, float(agent.get('E')))
    	co2 = np.append(co2, float(agent.get('c_CO2')))
    	co = np.append(co, float(agent.get('c_CO')))
    	hcl = np.append(hcl, float(agent.get('c_HCl')))
    	hcn = np.append(hcn, float(agent.get('c_HCN')))
    	fed_in = np.append(fed_in, float(agent.get('FED_In')))
    	fec_smoke = np.append(fec_smoke, float(agent.get('FEC_Smoke')))
    	fed_heat = np.append(fed_heat, float(agent.get('FED_Heat')))
    	fic_im = np.append(fic_im, float(agent.get('FIC_Im')))
    	fic_in = np.append(fic_in, float(agent.get('FED_In')))

### prepare for plotting
values = [E, co2, co, hcn, hcl, fed_in, fed_heat, fec_smoke, fic_im, fic_in]
labels = ["Extinction Coefficient (1/m)", "CO2 concentration (ppm)","CO concentration (ppm)",
 			"HCN concentration (ppm)", "HCl concentration (ppm)", "FED Incapacitation", "FED Heat Dose",
 			"FEC Smoke", "FIC Impairment", "FIC Incapacitation"]
units = ["1/m", "ppm", "ppm", "ppm", "ppm", "","","",""]

max_values = []

for i in values:
	max_values.append(max(i))
	i[i==0] = float('nan')

### extinction plot

plt.plot(t, E, 'bx')
plt.plot([0,max(t)], [max_values[0], max_values[0]], 'r-', label='max = %.2f %s' % (max_values[0], units[0]))
plt.xlabel("Time")
plt.ylabel(labels[0])
plt.ylim(0)
plt.legend()
plt.savefig('toxicity_output_plot_0_extinction.png')
plt.close()

### concentration plot
plt.figure(figsize=(12,9))

for i in range(1,5):
	plt.subplot(2,2,i)
	plt.plot(t, values[i], 'bx')
	plt.plot([0,max(t)], [max_values[i], max_values[i]], 'r-', label='max = %.2f %s' % (max_values[i], units[i]))
	plt.xlabel("Time")
	plt.ylabel(labels[i])
	plt.ylim(0)
	plt.legend()

plt.savefig('toxicity_output_plot_1_concentrations.png')
plt.close()

plt.figure(figsize=(12,9))
values2 = [fed_in, fec_smoke, fed_heat, fic_im, fic_in]

for j in range(5,9):
	plt.subplot(2,2,j-4)
	plt.plot(t, values[j], 'bx')
	plt.plot([0,max(t)], [max_values[j], max_values[j]], 'r-', label='max = %.3f %s' % (max_values[j], units[j]))
	plt.xlabel("Time")
	plt.ylabel(labels[j])
	plt.ylim(0)
	plt.legend()

plt.savefig('toxicity_output_plot_2_FED_FEC_FIC.png')
plt.close()