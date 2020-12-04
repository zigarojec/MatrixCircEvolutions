import matplotlib.pylab as plt
import runme2


a = runme2.evaluateActiveF_2_HP(9,21)
b = runme2.evaluateActiveF_2_HP(9,212)
c = runme2.evaluateActiveF_2_HP(9,2123)

plt.hold(True)
plt.title('Evolving filters Bode plot')
plt.ylabel('gain [dB]')
plt.xlabel('frequency [Hz]')

pa, = plt.semilogx(a['x']['nominal'], a['y']['nominal'], 'k--', linewidth=3, label='1st run')
pb, = plt.semilogx(b['x']['nominal'], b['y']['nominal'], 'k-', linewidth=3, label='2nd run (fine tuned)')
pc, = plt.semilogx(c['x']['nominal'], c['y']['nominal'], 'k-', linewidth=1, label='2nd run')

plt.legend(handles=[pa, pc, pb], loc=4)

plt.grid(True)
plt.hold(False)
plt.show()
