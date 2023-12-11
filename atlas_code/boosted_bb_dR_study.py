import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import ROOT

# Generate random numbers
n_events = 100000
theta_rand = np.random.uniform(0.0,3.1415,n_events)
#pH_rand = np.random.uniform(480, 7000,n_events)
pH_rand = []


# Code to generate falling exponential for pT of Higgs
pT_min = 480
pT_max = 7000

while len(pH_rand) < n_events:
  r = np.random.uniform(0, 1)
  x = -(1000) * np.log(1 - r)
  # Only keep x if it is in range [pT_min, pT_max]

  if ((x > pT_min) and (x < pT_max)):
    pH_rand.append(x)

# Constants
mH = 125.0
mB = 4.18
e1 = mH / 2.0
p1 = np.sqrt((np.power(mH,2) / 4.0) - np.power(mB,2))


counter = 0
delta_theta = []
fail_theta = []
fail_pH = []

for i,theta in enumerate(theta_rand):
  theta = theta_rand[i]
  pH = pH_rand[i]
  beta = pH / np.sqrt(np.power(mH,2) + np.power(pH,2))
  gamma = np.sqrt(np.power(mH,2) + np.power(pH,2)) / mH

  theta_prime_1 = np.arctan( (p1 * np.sin(theta)) / (gamma * ( p1 * np.cos(theta) - beta * e1)))
  theta_prime_2 = np.arctan( (p1 * np.sin(theta)) / (gamma * ( p1 * np.cos(theta) + beta * e1)))

#  theta_prime_1 = np.arctan((np.sqrt(np.power(mH,2)/2.0 - 3*np.power(mB,2)) * np.sin(theta)) / ((np.sqrt(np.power(mH,2) + np.power(pH,2)) / mH) * (np.sqrt(np.power(mH,2) - 3 * np.power( mB,2)) * np.cos(theta) - np.sqrt(np.power(mH,2) - np.power(mB,2)))))
#  theta_prime_2 = np.arctan( - (np.sqrt(np.power(mH,2)/2.0 - 3*np.power(mB,2)) * np.sin(theta)) / (((np.sqrt(np.power(mH,2) + np.power(pH,2))) / mH) * (np.sqrt(np.power(mH,2) - 3 * np.power( mB,2)) * np.cos(theta) + np.sqrt(np.power(mH,2) - np.power(mB,2)))))

  delta_theta.append(np.abs(theta_prime_1 - theta_prime_2))

  if(np.abs(theta_prime_1 - theta_prime_2) > 1.0):
    #print "Event with b falling outside R = 1.0 cone"
    #print "  theta = " + str(theta)
    #print "  higgs momentum = " + str(pH)
    fail_theta.append(theta)
    fail_pH.append(pH)
    counter += 1

print "With R = 1.0 cut we lose " + str(counter) + " out of " + str(n_events)

# Make figure with subplots
fig, axis_array = plt.subplots(2,2)

fig.suptitle("Study R = 1.0 fatjet ability to capture boosted b-jet daughters of higgs", fontsize = 16)

# the histogram of the data. hist(data, n_bins, normalization, color)
axis_array[0,0].hist(delta_theta, 100, color='green')
axis_array[0,0].set_xlabel('Phi, opening angle in Lab Frame [rad]')
axis_array[0,0].set_title("All Events (Green)")
axis_array[0,0].axis([0, 2.0, 0, 4000])
axis_array[0,0].grid(True)

axis_array[1,0].hist(pH_rand, 100, color='green')
axis_array[1,0].set_yscale('log', nonposy='clip')
axis_array[1,0].set_xlabel('pT of Higgs in Lab Frame [GeV/c]')
#axis_array[1,0].set_title("All events (Log)")
axis_array[1,0].axis([0, 7000, 0, 8000])
axis_array[1,0].grid(True)

axis_array[0,1].hist(fail_theta, 100, color='red')
#axis_array[1,0].hist(fail_theta, 25, color='red')
axis_array[0,1].set_xlabel('Theta, angle between Higgs and b-jets in Rest Frame [rad]')
axis_array[0,1].set_title("Events with R > 1.0 (Red)")
axis_array[0,1].axis([0, 3.1415, 0, 1000])
axis_array[0,1].grid(True)

axis_array[1,1].hist(fail_pH, 50, color='red')
axis_array[1,1].set_yscale('log', nonposy='clip')
axis_array[1,1].set_xlabel('Higgs momentum in Lab Frame [GeV/c]')
#axis_array[1,1].set_title("Evenets outside R = 1.0")
axis_array[1,1].axis([0, 3000, 0, 4000])
axis_array[1,1].grid(True)

#l = plt.plot(bins, y, 'r--', linewidth=1)

#plt.yscale('log')
#plt.xlabel('X-axis')
#plt.ylabel('Histogram')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.title("My Very First MatPlotLib Histogram")
#plt.axis([0, 2.0, 0, 1000])
#plt.grid(True)

plt.show()
