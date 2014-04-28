# -*- coding: latin1 -*-
#  Copyright (C) 2006 João Luís Silva <jsilva@fc.up.pt>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2, or (at your option)
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#

from __future__ import division

import gtk

try:
    import json
except ImportError:
    import simplejson as json

import scipy.interpolate as interpolate
import numpy as np
import numpy.fft as fft
from common import *
from widgets import *
import units as units_module
import glb

from graph_window import *
import os
#-----------------------------------------------------------------------
# Pulse shape
#-----------------------------------------------------------------------
class PulseShape:
    SHAPE_SECH,SHAPE_GAUSSIAN,SHAPE_SPECTRA = range(3)
    def __init__(self):
        self.pulse_shape = PulseShape.SHAPE_SECH
        self.graph = GraphWindow("x","y")
#        self.graph2 = GraphWindow("x","y")
        self.has_data = False   #Has spectra data, self.t,self.At,self.l,self.S

    def set_as_sech(self,t_fwhm):
        self.pulse_shape = PulseShape.SHAPE_SECH
        self.t_fwhm = t_fwhm
        self.t0 = t_fwhm / 1.7627

    def set_as_gaussian(self,t_fwhm):
        self.pulse_shape = PulseShape.SHAPE_GAUSSIAN
        self.t_fwhm = t_fwhm
        self.t0 = t_fwhm / 1.1774

    def set_as_spectra(self,l,S):
        self.pulse_shape = PulseShape.SHAPE_SPECTRA
        
        self.l = l
        self.S = S
        self.has_data = True

        #To calculate t0...
        l = np.array(l)*units.ToProgUnits(1.0,"nm")
        n = len(l)
        f0 = units.c/l[-1]
        f1 = units.c/l[0]
        dt = 1/(f1-f0)

        df = units.c/l[0]-units.c/l[1]

        nf = 2*n
        f = np.linspace(f0,f1,nf)
        Sf = np.linspace(f0,f1,nf)

        dt = 1/(f1-f0)
        df = f[1]-f[0]

        if dt>units.ToProgUnits(0.5,"fs"):
            target_F = 1/0.5
            n_new_points = (target_F - (f1-f0))/df
            Sf = np.concatenate((np.zeros(n_new_points/2),Sf))
            Sf = np.concatenate((Sf,np.zeros(n_new_points/2)))
            nf = len(Sf)
            f = np.linspace(f0-target_F/2.0,f0+target_F/2.0,nf)

        interp = interpolate.interp1d(l,S,bounds_error=False,fill_value=0.0)
        for i in range(nf):
            lambda_ = units.c/f[i]
            Sf[i] = interp(lambda_)*lambda_**2

        Sf = np.where(Sf<0.0,0.0,Sf)

#        print "CENTRAL FREQ",np.sum(Sf*f)/np.sum(Sf)

#        self.graph.ShowList(Sf,f,"f,Sf",None)
        At = fft.fftshift(fft.fft(np.sqrt(Sf)))
        It = np.abs(At)**2

        n = len(f)
        df = f[1]-f[0]
        F = n*df
        dt = 1/F
        T = n*dt
        t = np.linspace(-T/2.0,T/2.0,nf)

        self.t_fwhm = (t[1]-t[0])*CalculateFWHM(It)
        self.t0 =  self.t_fwhm / 1.1774

#        print "set_as_spectra: t0",self.t0


    def init(self,t,f):
        if self.pulse_shape != PulseShape.SHAPE_SPECTRA:
            return

        self.t = t
        self.freqs = f

        l = np.array(self.l)*units.ToProgUnits(1.0,"nm")
        m = np.amax(self.S)
        if m==0.0:
            m = 1.0
        S = self.S/m

        n = len(l)
        nf = len(f)

        #Convert to frequencies (f)
        Sf = np.zeros(nf)
        interp = interpolate.interp1d(l,S,bounds_error=False,fill_value=0.0)
        for i in range(nf):
            lambda_ = units.c/f[i]
        #    print "i",i,"f[i]",f[i],"lambda",lambda_
            Sf[i] = interp(lambda_)*lambda_**2

#        self.graph.ShowList(Sf,f,"f,Sf",None)
        Sf = np.where(Sf<0.0,0.0,Sf)

#        print "CENTRAL FREQ",np.sum(Sf*f)/np.sum(Sf)

#        self.graph.ShowList(Sf,f,"f,Sf",None)
        self.At = IFFT_t(np.sqrt(Sf))

#        self.At = IFFT_t(np.sqrt(Sf))
        It = np.abs(self.At)**2
        m = np.amax(It)
        if m!=0.0:
            It[:] /= m

        m = np.amax(np.abs(self.At))
        if m!=0.0:
            self.At[:] /= m

        self.t_fwhm = (self.t[1]-self.t[0])*CalculateFWHM(It)
        self.t0 =  self.t_fwhm / 1.1774


    def f(self,t,i):
        if self.pulse_shape == PulseShape.SHAPE_SECH:
            return 1.0/np.cosh(t/self.t0)
        elif self.pulse_shape == PulseShape.SHAPE_GAUSSIAN:
            return np.exp(-(t/self.t0)**2)
        elif self.pulse_shape == PulseShape.SHAPE_SPECTRA:
            if self.has_data:
                return self.At[i]

#-----------------------------------------------------------------------
# Simulation class
#-----------------------------------------------------------------------
class Simulation:
    TYPE_SELLMEIER,TYPE_TAYLOR,TYPE_D,TYPE_NONE=range(4)    #Dispersion Types
    TYPE_SOLITON,TYPE_GAUSSIAN,TYPE_FROM_SPECTRA,TYPE_P0,TYPE_PULSE_ENERGY=range(5)                       #Pulse types
    saved_dispersion_types = ["sellmeier","taylor","d","none"]
    saved_pulse_types = ["soliton","gaussian","spectra","p0","pulse_energy"]
    def __init__(self,prefs_):
        global prefs
        prefs = prefs_
        self.set_defaults()
        self.shape = PulseShape()

    def set_defaults(self):
        self.description = ""
        self.z_max = 0.0
        self.n_log = 50.0
        self.alpha = 0.0
        self.gamma = 0.0
        self.beta2 = 0.0
        self.disable_self_steepening = False
        self.type = Simulation.TYPE_NONE
        self.pulse_type = Simulation.TYPE_SOLITON
        self.pulse_int_type = Simulation.TYPE_SOLITON
        self.t_fwhm = 0.0
        self.gdd = 0.0
        self.tod = 0.0
        self.chirp2 = 0.0
        self.chirp3 = 0.0
        self.disable_Raman = False
        self.f_R = 0.18
        self.tau_1 = 12.2
        self.tau_2 = 32.0

        self.l = np.zeros(0,dtype=np.float64)
        self.D = np.zeros(0,dtype=np.float64)

        self.beta_center = 0.0
        self.beta = np.zeros(0,dtype=np.float64)

        self.ref_wavelength = 800.0
        self.w0 = l2w(units.ToProgUnits(self.ref_wavelength,"nm"))

        self.Sellmeier_B = np.zeros(0,dtype=np.float64)
        self.Sellmeier_l = np.zeros(0,dtype=np.float64)

    def Load(self,shortname):
        filename = os.path.join(prefs.data_dir,shortname+".sim")

        fp = open(filename,"r")
        try:
            m = json.load(fp)
        except:
            Error(glb.main.window,"An error occurred when opening the simulation "+shortname+" (file "+filename+")")
            return False

        fp.close()

        sim = Dict(m["simulation"])
        mat = Dict(m["material"])
        disp = Dict(m["dispersion"])

        self.set_defaults()

        self.description = sim.get("description","")
        self.z_max = sim.get("z_max",0.0)
        self.n_log = sim.get("n_log",50)

        if int(sim.get("disable_self_steepening",0)):
            self.disable_self_steepening = True
        else:
            self.disable_self_steepening = False

        if int(sim.get("disable_raman",0)):
            self.disable_Raman = True
            self.f_R = sim.get("f_r",0.18)
            self.tau_1 = sim.get("tau_1",12.2)
            self.tau_2 = sim.get("tau_1",32.0)
        else:
            self.disable_Raman = False

        strtype = sim.get("pulse_type","")
        try:
            self.pulse_type = Simulation.saved_pulse_types.index(strtype)
        except:
            Error(glb.simulation_editor.window,"Unknown pulse type")
            return False

        if self.pulse_type == Simulation.TYPE_SOLITON:
            self.t_fwhm = sim.get("t_fwhm",0.0)
            self.shape.set_as_sech(self.t_fwhm)
        elif self.pulse_type == Simulation.TYPE_GAUSSIAN:
            self.t_fwhm = sim.get("t_fwhm",0.0)
            self.shape.set_as_gaussian(self.t_fwhm)
        elif self.pulse_type == Simulation.TYPE_FROM_SPECTRA:
            self.spectra_l = sim["spectra_l"]
            self.spectra_S = sim["spectra_S"]
            self.shape.set_as_spectra(self.spectra_l,self.spectra_S)

        self.ref_wavelength = sim.get("ref_wavelength",800.0)
        self.w0 = l2w(units.ToProgUnits(self.ref_wavelength,"nm"))

        strtype = sim.get("pulse_int_type","")
        try:
            self.pulse_int_type = Simulation.saved_pulse_types.index(strtype)
        except:
            Error(glb.simulation_editor.window,"Unknown pulse intensity type")
            return False
        if self.pulse_int_type == Simulation.TYPE_SOLITON:
            self.soliton_order = sim.get("soliton_order",1.0)
        elif self.pulse_int_type == Simulation.TYPE_P0:
            self.p0 = sim.get("p0",1.0)
        else:
            self.pulse_energy = sim.get("pulse_energy",1.0)

        self.gdd = sim.get("gdd",0.0)
        self.tod = sim.get("tod",0.0)
        self.chirp2 = sim.get("chirp2",0.0)
        self.chirp3 = sim.get("chirp3",0.0)

        self.alpha = mat.get("alpha",0.0)
        self.gamma = mat.get("gamma",0.0)

        strtype = disp["type"]
        try:
            self.type = Simulation.saved_dispersion_types.index(strtype)
        except:
            Error(glb.simulation_editor.window,"Unknown dispersion type")
            return False


        if self.type == Simulation.TYPE_SELLMEIER:
            l = disp["sellmeier"]
            n = len(l)
            self.Sellmeier_B = np.zeros(n,dtype=np.float64)
            self.Sellmeier_l = np.zeros(n,dtype=np.float64)
            i = 0
            while i<n:
                self.Sellmeier_B[i] = l[i][0]
                self.Sellmeier_l[i] = l[i][1]
                i = i + 1

        elif self.type == Simulation.TYPE_TAYLOR:

            self.beta_center = disp["beta_center"]
            l = disp["betas"]
            n = len(l)
            self.beta = np.zeros(n,dtype=np.float64)

            i = 0
            while i<n:
                self.beta[i] = l[i]
                i = i + 1

        elif self.type == Simulation.TYPE_D:
            l = disp["d"]
            n = len(l)
            self.l = np.zeros(n,dtype=np.float64)
            self.D = np.zeros(n,dtype=np.float64)

            i = 0
            while i<n:
                self.l[i] = l[i][0]
                self.D[i] = l[i][1]
                i = i + 1
        elif self.type == Simulation.TYPE_NONE:
            self.beta2 = 0.0

        return True

    def CalculateBeta2(self,glb_omega):
        """Calculates beta(omega)-beta0 - beta1 * ( omega - omega0) where
        beta0 = beta(omega0) and beta1 = d beta / d omega at omega = omega0."""
        if self.type == Simulation.TYPE_NONE:
            self.betas = np.zeros(prefs.n,dtype=np.float64)
        elif self.type == Simulation.TYPE_D:
            #FIXME: units
            D = units.ToProgUnits(self.D,"ps nm^-1 km^-1")
            l = units.ToProgUnits(self.l,"um")

            #beta2 = -2 Pi c D(w0) / w0**2
            s = interpolate.InterpolatedUnivariateSpline(l,D)
            self.beta2 = -2.0*math.pi*units.c*s(w2l(self.w0))/self.w0**2

            betas = -D*(l**2)/(2.0*math.pi*units.c)
            
            omega = np.array([l2w(x) for x in l])
            omega.sort()
            betas = betas[::-1] #Revert

            s = interpolate.InterpolatedUnivariateSpline(omega,betas)

            omega2 = np.linspace(omega[0],omega[-1],len(omega)*20)
            n = len(omega2)
            I1 = np.zeros(n,dtype=np.float64)
            I2 = np.zeros(n,dtype=np.float64)

            w0 = self.w0

            for i in range(n):
                I1[i] = s.integral(w0,omega2[i])
#                print "Int(w0=",w0," to ",omega2[i],")=",I1[i]

            s = interpolate.InterpolatedUnivariateSpline(omega2,I1)

            for i in range(n):
                I2[i] = s.integral(w0,omega2[i])

            s = interpolate.InterpolatedUnivariateSpline(omega2,I2)
#            print omega2
            self.betas = s(glb_omega+w0)
        elif self.type == Simulation.TYPE_TAYLOR:
            self.betas = np.zeros(prefs.n,dtype=np.float64)
            w_c = 2.0*math.pi*units.c/units.ToProgUnits(self.beta_center,"nm")-self.w0
            omega_exp = glb_omega - w_c

            i = 2
            self.beta2 =  units.ToProgUnits(self.beta[0],"ps^2 km^-1")
            for b in self.beta:
                b_n = units.ToProgUnits(b,"ps^%d km^-1"%i)
                self.betas += b_n * (omega_exp**i) / Factorial(i)
                printf("CalculateBeta2, b = %g, b_n = %g,beta2 %g, betas %s\n",b,b_n,self.beta2,str(self.betas))
                i += 1
        else:
            raise LaserFOAM_Exception("CalculateBeta2: Not yet implemented!")

#-----------------------------------------------------------------------
# Simulation Editor
#-----------------------------------------------------------------------

class SimulationEditor:
    def __init__(self):
        global prefs
        Widgets(os.path.join("ui","simulation.ui")).connect(self)
        prefs = glb.prefs
        self.window = self.winSimulationEditor
        self.window.set_transient_for(glb.simulation_list.window)
        self.window.connect("delete_event",self.delete_event)

        self.winInt = GraphWindow("Time (fs)","Normalized intensity")
        self.winSpec = GraphWindow("Wavelength (nm)","Spectra")
        self.winD = GraphWindow("Lambda (micron)","Dispersion D (ps/(nm km))")

        self.response = gtk.RESPONSE_CANCEL

        self.ntbDispersion.set_show_tabs(False)
        self.ntbPulse.set_show_tabs(False)

        self.selection = self.tvwBetas.get_selection()
        self.selection.set_mode(gtk.SELECTION_BROWSE)
        
        self.liststore = gtk.ListStore(str,str,str)
        self.tvwBetas.set_model(self.liststore)

        self.cellCoeff = gtk.CellRendererText()
        self.cellValue = gtk.CellRendererText()
        self.cellUnits = gtk.CellRendererText()

        self.tvcolCoeff = gtk.TreeViewColumn('Coefficient',self.cellCoeff,text=0)
        self.tvcolValue = gtk.TreeViewColumn('Value',self.cellValue,text=1)
        self.tvcolUnits = gtk.TreeViewColumn('Units',self.cellUnits,text=2)

        self.tvwBetas.append_column(self.tvcolCoeff)
        self.tvwBetas.append_column(self.tvcolValue)
        self.tvwBetas.append_column(self.tvcolUnits)

        self.cellValue.set_property('editable',True)
        self.cellValue.connect('edited', self.cellValue_edited)
        self.imgInputPulse.set_from_file(os.path.join("images","input_pulse.png"))
        self.imgSellmeier.set_from_file(os.path.join("images","sellmeier.png"))
        self.imgSech.set_from_file(os.path.join("images","sech.png"))
        self.imgGaussian.set_from_file(os.path.join("images","gaussian.png"))
        self.imgRaman.set_from_file(os.path.join("images","raman.png"))

        self.selection_S = self.tvwSellmeier.get_selection()
        self.selection_S.set_mode(gtk.SELECTION_BROWSE)
        
        self.liststore_S = gtk.ListStore(str,str,str)
        self.tvwSellmeier.set_model(self.liststore_S)

        self.cellCoeff_S = gtk.CellRendererText()
        self.cellB = gtk.CellRendererText()
        self.cellLambda = gtk.CellRendererText()

        self.tvcolCoeff_S = gtk.TreeViewColumn('Coefficient',self.cellCoeff_S,text=0)
        self.tvcolB = gtk.TreeViewColumn('    B    ',self.cellB,text=1)
        self.tvcolLambda = gtk.TreeViewColumn('Lambda (micron)',self.cellLambda,text=2)

        self.tvwSellmeier.append_column(self.tvcolCoeff_S)
        self.tvwSellmeier.append_column(self.tvcolB)
        self.tvwSellmeier.append_column(self.tvcolLambda)

        self.cellB.set_property('editable',True)
        self.cellB.connect('edited', self.cellB_edited)

        self.cellLambda.set_property('editable',True)
        self.cellLambda.connect('edited', self.cellLambda_edited)

        self.selection_D = self.tvwD.get_selection()
        self.selection_D.set_mode(gtk.SELECTION_BROWSE)
        
        self.liststore_D = gtk.ListStore(float,float)
        self.tvwD.set_model(self.liststore_D)

        self.cellLambda_D = gtk.CellRendererText()
        self.cellValue_D = gtk.CellRendererText()

        self.tvcolLambda_D = gtk.TreeViewColumn('Lambda (micron)',self.cellLambda_D,text=0)
        self.tvcolValue_D = gtk.TreeViewColumn('D(lambda) (ps/(nm km))',self.cellValue_D,text=1)

        self.tvwD.append_column(self.tvcolLambda_D)
        self.tvwD.append_column(self.tvcolValue_D)

        self.sim = Simulation(prefs)
        self.cleanup()

        self.UpdateWindow()


    def delete_event(self,*args):
        """Do not destroy this window, just hide it"""
        self.window.set_property("visible",False)
        return True

    def on_simulation_editor_activate(self,*args):
        self.window.set_property("visible",True)
#        self.cleanup()

    def ProcessDispersion(self):
        if self.radioSellmeier.get_active():
            self.sim.type = Simulation.TYPE_SELLMEIER

            n = ListStoreSize(self.liststore_S)
            if n==0:
                Error(self.window,"The Sellmeier coefficients list is empty")
                return False

            B = []
            l = []
            for modelrow in self.liststore_S:
                v1 = self.liststore_S.get(modelrow.iter,1)[0]
                v2 = self.liststore_S.get(modelrow.iter,2)[0]
                if (not isFloat(v1)) or (not isFloat(v2)):
                    Error(self.window,"The Sellmeier coefficients must be numbers")
                    return False

                B.append(float(v1))
                l.append(float(v2))

            self.sim.Sellmeier_B = np.array(B,copy=True)
            self.sim.Sellmeier_l = np.array(l,copy=True)

        elif self.radioTaylor.get_active():
            self.sim.type = Simulation.TYPE_TAYLOR

            if not isFloat(self.txtMatBetaCenter.get_text()):
                Error(self.window,"The betas expansion center is invalid")
                return False

            self.sim.beta_center = float(self.txtMatBetaCenter.get_text())

            n = ListStoreSize(self.liststore)
            if n==0:
                Error(self.window,"The Taylor coefficients list is empty")
                return False

            b = []
            for modelrow in self.liststore:
                value = self.liststore.get(modelrow.iter,1)[0]
                if not isFloat(value):
                    Error(self.window,"The Taylor coefficients must be numbers")
                    return False

                b.append(float(value))

            self.sim.beta = np.array(b,copy=True)

        elif self.radioDisp.get_active():
            self.sim.type = Simulation.TYPE_D
            n = ListStoreSize(self.liststore_D)
            if n==0:
                Error(self.window,"The list is empty")
                return False

            l = []
            D = []
            for modelrow in self.liststore_D:
                l.append(self.liststore_D.get(modelrow.iter,0)[0])
                D.append(self.liststore_D.get(modelrow.iter,1)[0])

            self.sim.l = np.array(l,copy=True)
            self.sim.D = np.array(D,copy=True)
        else:
            self.sim.type = Simulation.TYPE_NONE

        self.UpdateWindow()
        return True

    def on_cmdSimulationEditorOk_clicked(self,*args):

        if not self.ProcessDispersion():
            return

        if not self.ValidateSimulation():
            return

        self.WriteSimulation(self.filename)
        glb.simulation_list.Refresh()
        self.window.set_property("visible",False)

#        self.DataToSim(self.sim)

#        glb.main.run(self.sim)


    def on_cmdSimulationEditorCancel_clicked(self,*args):
        self.window.set_property("visible",False)


    def ValidateSimulation(self):
        """Ensure the fields have valid values"""

        if not isFloat(self.txtZ_max.get_text()):
            Error(self.window,"The propagation length is not valid.")
            return False

        #Alpha and Gamma must be floating point values
        if not isFloat(self.txtMatAlpha.get_text()):
            Error(self.window,"Alpha is not valid.")
            return False

        if not isFloat(self.txtMatGamma.get_text()):
            Error(self.window,"Gamma is not valid.")
            return False

        if not isFloat(self.txtRefWavelength.get_text()):
            Error(self.window,"The reference wavelength is not valid.")
            return False

        if float(self.txtRefWavelength.get_text())<=0.0:
            Error(self.window,"The reference wavelength is not valid.")
            return False

        if not isFloat(self.txtF_R.get_text()):
            Error(self.window,"f_R is not valid.")
            return False

        if float(self.txtF_R.get_text())<=0.0 or float(self.txtF_R.get_text())>1.0:
            Error(self.window,"f_R must be between 0.0 and 1.0")
            return False

        if not isFloat(self.txtTau1.get_text()):
            Error(self.window,"Tau1 is not valid.")
            return False

        if not isFloat(self.txtTau2.get_text()):
            Error(self.window,"Tau2 is not valid.")
            return False

        if self.radioSolitonOrder.get_active():
            if self.radioNoDisp.get_active():
                Error(self.window,"You cannot specify the soliton order when there is no dispersion")
                return False

            if float(self.txtMatGamma.get_text())==0.0:
                Error(self.window,"You cannot specify the soliton order when the nonlinear coefficient gamma is zero")
                return False

            if not isFloat(self.txtSolitonOrder.get_text()):
                Error(self.window,"The soliton order is not valid.")
                return False

            if float(self.txtSolitonOrder.get_text())<=0.0:
                Error(self.window,"The soliton order is not valid.")
                return False
        elif self.radioPeakPower.get_active():
            if not isFloat(self.txtPeakPower.get_text()):
                Error(self.window,"The peak power is not valid.")
                return False

            if float(self.txtPeakPower.get_text())<=0.0:
                Error(self.window,"The peak power is not valid.")
                return False
        else:
            if not isFloat(self.txtPulseEnergy.get_text()):
                Error(self.window,"The pulse energy is not valid.")
                return False

            if float(self.txtPulseEnergy.get_text())<=0.0:
                Error(self.window,"The pulse energy is not valid.")
                return False


        if not isFloat(self.txtGDD.get_text()):
            Error(self.window,"The GDD is not valid.")
            return False

        if not isFloat(self.txtTOD.get_text()):
            Error(self.window,"The TOD is not valid.")
            return False

        if not isFloat(self.txtChirp2.get_text()):
            Error(self.window,"The chirp parameter C is not valid.")
            return False

        if not isFloat(self.txtChirp3.get_text()):
            Error(self.window,"The chirp parameter D is not valid.")
            return False

        if self.radioPulseSoliton.get_active():
            if not isFloat(self.txtSolitonDuration.get_text()):
                Error(self.window,"The duration is not valid.")
                return False

            if float(self.txtSolitonDuration.get_text()) <=0.0:
                Error(self.window,"The duration is not valid.")
                return False

        elif self.radioPulseGaussian.get_active():
            if not isFloat(self.txtGaussianDuration.get_text()):
                Error(self.window,"The duration is not valid.")
                return False

            if float(self.txtGaussianDuration.get_text()) <=0.0:
                Error(self.window,"The duration is not valid.")
                return False

        elif self.radioPulseSpectra.get_active():
            if not self.sim.shape.has_data:
                Error(self.window,"You haven't specified the spectra to import yet.")
                return False

        return True

    def WriteSimulation(self,filename):
        """Creates a simulation file"""
        self.DataToSim(self.sim)

        mat = {} 
        disp = {}

        mat["alpha"] = float(self.txtMatAlpha.get_text())
        mat["gamma"] = float(self.txtMatGamma.get_text())

        disp["type"] = Simulation.saved_dispersion_types[self.sim.type]


        if self.sim.type == Simulation.TYPE_SELLMEIER:
            disp["sellmeier"] = [[self.sim.Sellmeier_B[j],self.sim.Sellmeier_l[j]] for j in range(len(self.sim.Sellmeier_B))]
        elif self.sim.type == Simulation.TYPE_TAYLOR:
            disp["beta_center"] = self.sim.beta_center
            disp["betas"]=[e for e in self.sim.beta]
        elif self.sim.type == Simulation.TYPE_D:
            n = len(self.sim.l)
            disp["d"] = [[self.sim.l[j],self.sim.D[j]] for j in range(n)]

        sim = {}
        sim["description"] = self.sim.description
        sim["z_max"] = self.sim.z_max
        sim["n_log"] = self.sim.n_log
        sim["disable_raman"] = int(self.sim.disable_Raman)
        if not int(self.sim.disable_Raman):
            sim["f_r"] = self.sim.f_R
            sim["tau_1"] = self.sim.tau_1
            sim["tau_2"] = self.sim.tau_2
        sim["pulse_type"] = Simulation.saved_pulse_types[self.sim.pulse_type]
        sim["disable_self_steepening"] = int(self.sim.disable_self_steepening)
        sim["ref_wavelength"] = self.sim.ref_wavelength

        sim["pulse_int_type"] = Simulation.saved_pulse_types[self.sim.pulse_int_type]
        if self.sim.pulse_int_type == Simulation.TYPE_SOLITON:
            sim["soliton_order"] = self.sim.soliton_order
        elif self.sim.pulse_int_type == Simulation.TYPE_P0:
            sim["p0"] = self.sim.p0
        else:
            sim["pulse_energy"] = self.sim.pulse_energy

        if self.sim.pulse_type in (Simulation.TYPE_SOLITON,Simulation.TYPE_GAUSSIAN):
            sim["t_fwhm"] = self.sim.t_fwhm
        elif self.sim.pulse_type == Simulation.TYPE_FROM_SPECTRA:
            sim["spectra_l"] = list(self.sim.shape.l)
            sim["spectra_S"] = list(self.sim.shape.S[:])


        sim["gdd"] = self.sim.gdd
        sim["tod"] = self.sim.tod
        sim["chirp2"] = self.sim.chirp2
        sim["chirp3"] = self.sim.chirp3

        fp = open(filename,"w")
        json.dump({"simulation":sim,"material":mat,"dispersion":disp},fp,indent=1)
        fp.close()

    def New(self,shortname):
        """Creates a new simulation"""
        self.filename = os.path.join(prefs.data_dir,shortname+".sim")
        self.cleanup()
        self.window.set_property("visible",True)

    def Edit(self,shortname):
        """Edit the specified simulation"""
        self.filename = os.path.join(prefs.data_dir,shortname+".sim")
        self.cleanup()
        if not self.ReadSimulation(shortname):
            return
        self.window.set_property("visible",True)

    def ReadDispersion(self,sim):
        if sim.type == Simulation.TYPE_SELLMEIER:
            self.radioSellmeier.set_active(True)
            self.liststore_S.clear()
            n = len(self.sim.Sellmeier_B)
            i = 0
            while i<n:
                self.liststore_S.append((i+1,sim.Sellmeier_B[i],sim.Sellmeier_l[i]))
                i = i + 1
        elif sim.type == Simulation.TYPE_TAYLOR:
            self.radioTaylor.set_active(True)
            self.txtMatBetaCenter.set_text(str(sim.beta_center))

            n = len(sim.beta)
            i = 0
            while i<n:
                self.liststore.append((str(i+2),sim.beta[i],"ps^%d/km" % (i+2)))
                i = i + 1
        elif sim.type == Simulation.TYPE_D:
            self.radioDisp.set_active(True)
            n = len(sim.D)
            i = 0
            while i<n:
                self.liststore_D.append((sim.l[i],sim.D[i]))
                i = i + 1

    def ReadSimulation(self,shortname):
        """Reads a simulation file"""
        self.sim = Simulation(prefs)
        sim = self.sim
        if not sim.Load(shortname):
            return False

        self.txtSimulationDesc.set_text(sim.description)
        self.txtZ_max.set_text(str(sim.z_max))
        self.txtN_Log.set_text(str(sim.n_log))


        if sim.pulse_type == Simulation.TYPE_SOLITON:
            self.radioPulseSoliton.set_active(True)
            self.txtSolitonDuration.set_text(str(sim.t_fwhm))
            self.txtGaussianDuration.set_text(str(sim.t_fwhm))
        elif sim.pulse_type == Simulation.TYPE_GAUSSIAN:
            self.radioPulseGaussian.set_active(True)
            self.txtSolitonDuration.set_text(str(sim.t_fwhm))
            self.txtGaussianDuration.set_text(str(sim.t_fwhm))
        elif sim.pulse_type == Simulation.TYPE_FROM_SPECTRA:
            self.radioPulseSpectra.set_active(True)

        if sim.pulse_int_type == Simulation.TYPE_SOLITON:
            self.txtSolitonOrder.set_text(str(sim.soliton_order))
            self.radioSolitonOrder.set_active(True)
        elif sim.pulse_int_type == Simulation.TYPE_P0:
            self.txtPeakPower.set_text(str(sim.p0))
            self.radioPeakPower.set_active(True)
        else:
            self.txtPulseEnergy.set_text(str(sim.pulse_energy))
            self.radioPulseEnergy.set_active(True)

        self.txtGDD.set_text(str(sim.gdd))
        self.txtTOD.set_text(str(sim.tod))
        self.txtChirp2.set_text(str(sim.chirp2))
        self.txtChirp3.set_text(str(sim.chirp3))

#        self.D = mat.D

        self.txtMatAlpha.set_text(str(sim.alpha))
        self.txtMatGamma.set_text(str(sim.gamma))
        self.txtRefWavelength.set_text(str(sim.ref_wavelength))

        self.chkDisableSelfSteepening.set_active(sim.disable_self_steepening)
        self.chkDisableRaman.set_active(sim.disable_Raman)
        if not sim.disable_Raman:
            self.txtF_R.set_text(str(sim.f_R))
            self.txtTau1.set_text(str(sim.tau_1))
            self.txtTau2.set_text(str(sim.tau_2))

        self.ReadDispersion(sim)

        return True

    def DataToSim(self,sim):
        sim.description = self.txtSimulationDesc.get_text()
        sim.z_max = float(self.txtZ_max.get_text())
        sim.n_log = float(self.txtN_Log.get_text())
        sim.f_R = float(self.txtF_R.get_text())
        sim.tau_1 = float(self.txtTau1.get_text())
        sim.tau_2 = float(self.txtTau2.get_text())
        sim.disable_self_steepening = self.chkDisableSelfSteepening.get_active()
        sim.disable_Raman = self.chkDisableRaman.get_active()
        sim.ref_wavelength = float(self.txtRefWavelength.get_text())

        if self.radioPulseSoliton.get_active():
            sim.pulse_type = Simulation.TYPE_SOLITON
            sim.t_fwhm = float(self.txtSolitonDuration.get_text())
        elif self.radioPulseGaussian.get_active():
            sim.pulse_type = Simulation.TYPE_GAUSSIAN
            sim.t_fwhm = float(self.txtGaussianDuration.get_text())
        elif self.radioPulseSpectra.get_active():
            sim.pulse_type = Simulation.TYPE_FROM_SPECTRA

        sim.gdd = float(self.txtGDD.get_text())
        sim.tod = float(self.txtTOD.get_text())
        sim.chirp2 = float(self.txtChirp2.get_text())
        sim.chirp3 = float(self.txtChirp3.get_text())
        if self.radioSolitonOrder.get_active():
            sim.pulse_int_type = Simulation.TYPE_SOLITON
            sim.soliton_order = float(self.txtSolitonOrder.get_text())
        elif self.radioPeakPower.get_active():
            sim.pulse_int_type = Simulation.TYPE_P0
            sim.p0 = float(self.txtPeakPower.get_text())
        else:
            sim.pulse_int_type = Simulation.TYPE_PULSE_ENERGY
            sim.pulse_energy = float(self.txtPulseEnergy.get_text())

    def cleanup(self):
        """Initialize all the fields to their default values"""
        self.txtZ_max.set_text("")
        self.txtN_Log.set_text("50")
        self.txtF_R.set_text("0.18")
        self.txtTau1.set_text("12.2")
        self.txtTau2.set_text("32.0")
        self.txtSimulationDesc.set_text("")
        self.txtMatAlpha.set_text("")
        self.txtMatGamma.set_text("")
        self.txtRefWavelength.set_text("")
        self.txtPeakPower.set_text("")
        self.txtPulseEnergy.set_text("")
        self.txtZ_max.set_text("")
        self.liststore.clear()
        self.liststore_S.clear()
        self.liststore_D.clear()
        self.txtMatBetaCenter.set_text("")
        self.txtSolitonOrder.set_text("")
        self.txtSolitonDuration.set_text("")
        self.txtGaussianDuration.set_text("")
        self.txtGDD.set_text("0")
        self.txtTOD.set_text("0")
        self.txtChirp2.set_text("0")
        self.txtChirp3.set_text("0")
        self.lblSpectralWidth.set_text("")
        self.radioPulseSoliton.set_active(True)
        self.chkDisableSelfSteepening.set_active(False)
        self.chkDisableRaman.set_active(False)
        self.ntbSimulationEditor.set_current_page(0)
        self.radioNoDisp.set_active(True)
        self.radioPulseSoliton.set_active(True)
        self.radioSolitonOrder.set_active(True)
        self.sim.shape.has_data = False

    def UpdateWindow(self):
        """Enables or disables the widgets according to the selected options"""
#        print "UpdateWindow: self.sim.shape.has_data?",self.sim.shape.has_data
        if self.radioSellmeier.get_active():
            self.ntbDispersion.set_current_page(0)
        elif self.radioTaylor.get_active():
            self.ntbDispersion.set_current_page(1)
        elif self.radioDisp.get_active():
            self.ntbDispersion.set_current_page(2)
        elif self.radioNoDisp.get_active():
            self.ntbDispersion.set_current_page(3)

        if self.radioSolitonOrder.get_active():
            self.txtSolitonOrder.set_sensitive(True)
            self.txtPeakPower.set_sensitive(False)
            self.txtPulseEnergy.set_sensitive(False)
        elif self.radioPeakPower.get_active():
            self.txtSolitonOrder.set_sensitive(False)
            self.txtPeakPower.set_sensitive(True)
            self.txtPulseEnergy.set_sensitive(False)
        else:
            self.txtSolitonOrder.set_sensitive(False)
            self.txtPeakPower.set_sensitive(False)
            self.txtPulseEnergy.set_sensitive(True)

        if self.radioPulseSoliton.get_active():
            self.ntbPulse.set_current_page(0)
        elif self.radioPulseGaussian.get_active():
            self.ntbPulse.set_current_page(1)
        elif self.radioPulseSpectra.get_active():
            self.ntbPulse.set_current_page(2)
            self.cmdShowIntensity.set_sensitive(self.sim.shape.has_data)
            self.cmdShowSpectra.set_sensitive(self.sim.shape.has_data)

        if self.chkDisableRaman.get_active():
            self.frmRaman.set_sensitive(False)
            if self.chkDisableSelfSteepening.get_active():
                self.imgEquation.set_from_file(os.path.join("images","nls_no_ss.png"))
            else:
                self.imgEquation.set_from_file(os.path.join("images","nls.png"))
        else:
            self.frmRaman.set_sensitive(True)
            if self.chkDisableSelfSteepening.get_active():
                self.imgEquation.set_from_file(os.path.join("images","nls_raman_no_ss.png"))
            else:
                self.imgEquation.set_from_file(os.path.join("images","nls_raman.png"))

    def DurationAndSpectralWidth(self,t0,disp2,disp3,chirp2,chirp3):
        """Returns the duration and spectral width (in frequency units) for a soliton with the given parameters."""
        n = 4096
        if t0<=0.0:
            return -1,-1
        A = np.zeros(n,dtype=np.complex128)
        T = math.pi*t0*math.sqrt(n)
        t = np.zeros(n,dtype=np.float64)
        dt = T / (n-1)
        for i in range(n):
            t[i] = -T/2.0+i*dt
        omega = self.sim.w0 + CalculateOmega(n,T)
        f = omega/(2.0*np.pi)

        self.sim.shape.init(t,f)

        for i in range(n):
            A[i] = self.sim.shape.f(t[i],i)*np.exp(-1j*(chirp2/2.0)*(t[i]/t0)**2-1j*(chirp3/3.0)*(t[i]/t0)**3)

        A = IFFT_t(np.exp(1j*(disp2/2.0)*omega**2+1j*(disp3/6.0)*omega**3)*FFT_t(A))

        t_fwhm = CalculateFWHM(np.abs(A)**2)*dt

        S = np.abs(FFT_t(A))**2
        fwhm = CalculateFWHM(S)
        if fwhm == -1:
            return -1,-1

        return t_fwhm,fwhm/T

    def on_chkDisableSelfSteepening_clicked(self,*args):
        self.UpdateWindow()

    def on_chkDisableRaman_clicked(self,*args):
        self.UpdateWindow()


    def cellValue_edited(self,cell_renderer,path,new_value):
        iter = self.liststore[int(path)].iter
        self.liststore.set_value(iter,1,str(float(new_value)))

    def cellB_edited(self,cell_renderer,path,new_value):
        iter = self.liststore_S[int(path)].iter
        self.liststore_S.set_value(iter,1,str(float(new_value)))

    def cellLambda_edited(self,cell_renderer,path,new_value):
        iter = self.liststore_S[int(path)].iter
        self.liststore_S.set_value(iter,2,str(float(new_value)))

    def on_dispersion_editor_radiobutton_toggled(self,*args):
        self.UpdateWindow()

    def on_radioPulseType_toggled(self,*args):
        self.UpdateWindow()
        self.update_lblSpectralWidth() 

    def on_radioPulseIntensity_toggled(self,*args):
        self.UpdateWindow()

    def on_cmdBetasAdd_clicked(self,*args):
        n = ListStoreSize(self.liststore)
        self.liststore.append((str(n+2),"0.0","ps^%d/km" % (n+2)))

    def on_cmdBetasRemove_clicked(self,*args):
        model,l = self.selection.get_selected_rows()
        if len(l)!=1:
            return                  #FIXME: No row selected - give some feedback ?

        if ListStoreSize(self.liststore)==1:
            Message(self.window,"Cannot remove the last coefficient")
            return

        ind = l[0][0]
        iter = self.liststore[ind].iter
        self.liststore.remove(iter)

        #Renumber the list
        i = 2
        for modelrow in self.liststore:
            iter = modelrow.iter
            self.liststore.set_value(iter,0,str(i))
            self.liststore.set_value(iter,2,"ps^%d/km" % i)
            i = i + 1

    def on_cmdSellmeierAdd_clicked(self,*args):
        n = ListStoreSize(self.liststore_S)
        self.liststore_S.append((str(n+1),"0.0","0.0"))

    def on_cmdSellmeierRemove_clicked(self,*args):
        model,l = self.selection_S.get_selected_rows()
        if len(l)!=1:
            return                  #FIXME: No row selected - give some feedback ?

        if ListStoreSize(self.liststore_S)==1:
            Message(self.window,"Cannot remove the last coefficient")
            return

        ind = l[0][0]
        iter = self.liststore_S[ind].iter
        self.liststore_S.remove(iter)

        #Renumber the list
        i = 1
        for modelrow in self.liststore_S:
            iter = modelrow.iter
            self.liststore_S.set_value(iter,0,str(i))
            i = i + 1

    def on_cmdImportMatDispersion_clicked(self,*args):
        dialog = gtk.FileChooserDialog("Open..",
                                   None,
                                   gtk.FILE_CHOOSER_ACTION_OPEN,
                                   (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                    gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        response = dialog.run()
        if response == gtk.RESPONSE_OK:
            filename = dialog.get_filename()
            if not os.path.exists(filename):
                Error(glb.main.window,"The file "+filename+" doesn't exist")
                dialog.destroy()
                return False

            data = LoadCSV(filename)
            self.liststore_D.clear()
            for e in data:
                self.liststore_D.append((e[0],e[1]))
            
        dialog.destroy()

    def on_cmdShowMatDispersion_clicked(self,*args):
        l = []
        D = []
        for modelrow in self.liststore_D:
            l.append(self.liststore_D.get(modelrow.iter,0)[0])
            D.append(self.liststore_D.get(modelrow.iter,1)[0])
        self.winD.ShowList(D,l,"Dispersion",None)

    def on_txtRefWavelength_changed(self,*args):
        self.update_lblSpectralWidth() 

    def on_txtSolitonDuration_changed(self,*args):
        self.update_lblSpectralWidth() 

    def on_txtGaussianDuration_changed(self,*args):
        self.update_lblSpectralWidth() 

    def on_txtChirpParameter_changed(self,*args):
        self.update_lblSpectralWidth() 

    def update_lblSpectralWidth(self,*args):
        if not isFloat(self.txtGDD.get_text()) or not isFloat(self.txtTOD.get_text()) or not isFloat(self.txtChirp2.get_text()) or not isFloat(self.txtChirp3.get_text()) or not isFloat(self.txtRefWavelength.get_text()):
            self.lblSpectralWidth.set_text("Spectral width: ?? nm ( ?? THz)")
            return

        if self.radioPulseSoliton.get_active():
            if not isFloat(self.txtSolitonDuration.get_text()):
                self.lblSpectralWidth.set_text("Spectral width: ?? nm ( ?? THz)")
                return
            t_fwhm = float(self.txtSolitonDuration.get_text())
            self.sim.shape.set_as_sech(t_fwhm)
        elif self.radioPulseGaussian.get_active():
            if not isFloat(self.txtGaussianDuration.get_text()):
                self.lblSpectralWidth.set_text("Spectral width: ?? nm ( ?? THz)")
                return
            t_fwhm = float(self.txtGaussianDuration.get_text())
            self.sim.shape.set_as_gaussian(t_fwhm)
        elif self.radioPulseSpectra.get_active():
            if not self.sim.shape.has_data:
                self.lblSpectralWidth.set_text("Spectral width: ?? nm ( ?? THz)")
                return

        t0 = self.sim.shape.t0
        GDD = float(self.txtGDD.get_text())
        TOD = float(self.txtTOD.get_text())
        chirp2 = float(self.txtChirp2.get_text())
        chirp3 = float(self.txtChirp3.get_text())
        ref_wavelength = units.ToProgUnits(float(self.txtRefWavelength.get_text()),"nm")

        if ref_wavelength <= 0.0:
            self.lblSpectralWidth.set_text("Spectral width: ?? nm ( ?? THz)")
            return

        self.lblSpectralWidth.set_text("t0 %g GDD %g TOD %g ref_wavelength %g" % (t0,GDD,TOD,ref_wavelength))
        duration,spectral_width = self.DurationAndSpectralWidth(t0,GDD,TOD,chirp2,chirp3)

        if spectral_width == -1:
            self.lblSpectralWidth.set_text("Spectral width: ?? nm ( ?? THz)")
            return

        w0 = l2w(ref_wavelength)
        dl = w2l(w0-2.0*math.pi*spectral_width/2.0)-w2l(w0+2.0*math.pi*spectral_width/2.0)

        self.lblSpectralWidth.set_text("Spectral width: %g nm ( %g THz) Duration (FWHM): %g fs" % (units.FromProgUnits(dl,"nm"),units.FromProgUnits(spectral_width,"s^-1")/1.0E12,units.FromProgUnits(duration,"fs")))

    def on_cmdImportSpectra_clicked(self,*args):
        dialog = gtk.FileChooserDialog("Open..",
                                   None,
                                   gtk.FILE_CHOOSER_ACTION_OPEN,
                                   (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                    gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        response = dialog.run()
        if response == gtk.RESPONSE_OK:
            filename = dialog.get_filename()
            if not os.path.exists(filename):
                Error(glb.main.window,"The file "+filename+" doesn't exist")
                dialog.destroy()
                return False
            l,S = LoadData(filename)
            self.sim.shape.set_as_spectra(l,S)
            self.UpdateWindow()
            self.update_lblSpectralWidth() 
        dialog.destroy()

    def on_cmdShowIntensity_clicked(self,*args):
        self.winInt.ShowList(np.abs(self.sim.shape.At)**2,self.sim.shape.t,"Calculated intensity from the imported spectra")

    def on_cmdShowSpectra_clicked(self,*args):
        self.winSpec.ShowList(self.sim.shape.S,self.sim.shape.l,"Imported spectra")

units = units_module.Units()
