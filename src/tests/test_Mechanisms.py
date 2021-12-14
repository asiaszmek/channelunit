import os
from subprocess import run
import unittest

from neuron import h
import neuron
from channelunit import mechanisms_path

working_dir = os.getcwd()
os.chdir(mechanisms_path)
p = run('nrnivmodl')
neuron.load_mechanisms(mechanisms_path)
os.chdir(working_dir)
h.load_file("stdrun.hoc")


class TestMyVoltageClamp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.soma = h.Section("soma")
        cls.clamp = h.SEClampOLS(cls.soma(0.5))
        cls.clamp.dur1 = 10
        cls.clamp.amp1 = 100
        cls.clamp.dur2 = 20
        cls.clamp.amp2 = 200
        cls.clamp.dur3 = 30
        cls.clamp.amp3 = 300
        cls.clamp.dur4 = 40
        cls.clamp.amp4 = 400
        cls.clamp.dur5 = 50
        cls.clamp.amp5 = 500
        cls.clamp.dur6 = 60
        cls.clamp.amp6 = 600
        cls.clamp.dur7 = 70
        cls.clamp.amp7 = 700
        cls.clamp.dur8 = 80
        cls.clamp.amp8 = 800
        cls.clamp.dur9 = 90
        cls.clamp.amp9 = 900
        cls.clamp.dur10 = 100
        cls.clamp.amp10 = 1000
        cls.clamp.dur11 = 110
        cls.clamp.amp11 = 1100
        cls.clamp.dur12 = 120
        cls.clamp.amp12 = 1200
        cls.dt = 1
        h.cvode_active(1)
        time = h.Vector()
        time.record(h._ref_t, cls.dt)
        cls.voltage = h.Vector()
        cls.voltage.record(cls.clamp._ref_vc, cls.dt)
        h.init()
        h.tstop = (10+120)*12/2+10
        h.run()


    def test_amp1(self):
        out = set(self.voltage.as_numpy()[:int(10/self.dt)].tolist())
        self.assertEqual(out, set([100.]))

    def test_amp2(self):
        out = set(self.voltage.as_numpy()[int(10/self.dt):int(30/self.dt)].tolist())
        self.assertEqual(out, set([200.]))

    def test_amp3(self):
        out = set(self.voltage.as_numpy()[int(30/self.dt):int(60/self.dt)].tolist())
        self.assertEqual(out, set([300.]))

    def test_amp4(self):
        out = set(self.voltage.as_numpy()[int(60/self.dt):int(100/self.dt)].tolist())
        self.assertEqual(out, set([400.]))

    def test_amp5(self):
        out = set(self.voltage.as_numpy()[int(100/self.dt):int(150/self.dt)].tolist())
        self.assertEqual(out, set([500.]))

    def test_amp6(self):
        out = set(self.voltage.as_numpy()[int(150/self.dt):int(210/self.dt)].tolist())
        self.assertEqual(out, set([600.]))

    def test_amp7(self):
        out = set(self.voltage.as_numpy()[int(210/self.dt):int(280/self.dt)].tolist())
        self.assertEqual(out, set([700.]))

    def test_amp8(self):
        out = set(self.voltage.as_numpy()[int(280/self.dt):int(360/self.dt)].tolist())
        self.assertEqual(out, set([800.]))

    def test_amp9(self):
        out = set(self.voltage.as_numpy()[int(360/self.dt):int(450/self.dt)].tolist())
        self.assertEqual(out, set([900.]))

    def test_amp10(self):
        out = set(self.voltage.as_numpy()[int(450/self.dt):int(550/self.dt)].tolist())
        self.assertEqual(out, set([1000.]))

    def test_amp11(self):
        out = set(self.voltage.as_numpy()[int(550/self.dt):int(660/self.dt)].tolist())
        self.assertEqual(out, set([1100.]))

    def test_amp12(self):
        out = set(self.voltage.as_numpy()[int(660/self.dt):int(780/self.dt)].tolist())
        self.assertEqual(out, set([1200.]))

if __name__ == "__main__":
    unittest.main()
