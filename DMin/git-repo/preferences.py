# -*- coding: latin1 -*-
#  Copyright (C) 2006-2009 João Luís Silva <jsilva@fc.up.pt>
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

from common import *

import ConfigParser
import codecs
import math
import os.path
import units as units_module

from widgets import *

#Imports specific to MS Windows 
import platform
if platform.uname()[0].lower()=="windows":
    from win32com.shell import shellcon, shell

def HomeDir():
    """Returns the user home directory. Works in Linux and Windows"""
    if platform.uname()[0].lower()=="windows":
        return shell.SHGetFolderPath(0, shellcon.CSIDL_PROFILE, 0, 0)
    else:
        return os.path.expanduser("~")


#-----------------------------------------------------------------------
# Preferences
#-----------------------------------------------------------------------

class Preferences:
    def __init__(self,home_dir=None,lf_dir=None,data_dir=None):
        """Initialize the preferences from the laserfoamrc file. The default paths can be overridden."""
        self.created_rc = False
        if home_dir:
            self.home_dir = home_dir
        else:
            self.home_dir = HomeDir()

        self.lf_dir = lf_dir
        self.data_dir = data_dir

        self.Load()

    def Load(self):
        """Load the preferences file laserfoamrc. If it does not exist, create a default one"""
        if not self.lf_dir:
            self.lf_dir = os.path.join(self.home_dir,".laserfoam")
        self.rc_filename = os.path.join(self.lf_dir,"laserfoamrc")

        if os.path.exists(self.rc_filename):
            config = ConfigParser.ConfigParser()
            fp = codecs.open(self.rc_filename,"rb",encoding="utf-8")
            config.readfp(fp)
            if not config.has_section("prefs"):
                Error(None,"The preferences file "+self.rc_filename +" must have a section named [prefs].")

            if not config.has_option("prefs","data_dir"):
                Error(None,"The preferences file "+self.rc_filename +" must have a key called data_dir on the section [prefs].")

            if not self.data_dir:
                self.data_dir = config.get("prefs","data_dir")
                if not os.path.exists(self.data_dir):
                    Error(None,"The data path \""+ self.data_dir +"\" specified on the preferences file \""+self.rc_filename +"\" does not exist.")
            else:
                if not os.path.exists(self.data_dir):
                    Error(None,"The specified data path \""+ self.data_dir +"\" does not exist.")

            if not config.has_option("prefs","maximize"):
                config.set("prefs","maximize","1")

            if int(config.get("prefs","maximize")):
                self.maximize = True
            else:
                self.maximize = False

            if not config.has_option("prefs","2n"):
                config.set("prefs","2n","12")

            self.power = int(config.get("prefs","2n"))
            self.n = 2**self.power

            if not config.has_option("prefs","local_error"):
                config.set("prefs","local_error","0.001")

            self.local_error = float(config.get("prefs","local_error"))

            if not config.has_option("prefs","update_interval"):
                config.set("prefs","update_interval","2.0")

            self.update_interval = float(config.get("prefs","update_interval"))
            if not config.has_option("prefs","graphics"):
                config.set("prefs","graphics","wavelength")

            if not config.has_option("prefs","lambda0"):
                config.set("prefs","lambda0","100")

            if not config.has_option("prefs","lambda1"):
                config.set("prefs","lambda1","2000")

            self.lambda0 = config.get("prefs","lambda0")
            self.lambda1 = config.get("prefs","lambda1")

            self.freq_min = units.c/units.ToProgUnits(float(self.lambda1),"nm")
            self.freq_max = units.c/units.ToProgUnits(float(self.lambda0),"nm")

            if str(config.get("prefs","graphics")).lower()=="freq":
                self.use_wavelength = False
            else:
                self.use_wavelength = True

            fp.close()
        else: #laserfoamrc does not exists. Create it.

            #Does ~/.laserfoam exists?
            if not os.path.exists(os.path.join(self.home_dir,".laserfoam")):
                os.mkdir(os.path.join(self.home_dir,".laserfoam"))

            self.data_dir = os.path.join(self.home_dir,".laserfoam","data")
            if not os.path.exists(self.data_dir):
                os.mkdir(self.data_dir)

            self.maximize = True
            self.update_interval = 2.0
            self.local_error = 0.0001
            self.power = 12
            self.n = 2**self.power
            self.lambda0 = 100
            self.lambda1 = 2000
            self.freq_min = units.c/units.ToProgUnits(float(self.lambda1),"nm")
            self.freq_max = units.c/units.ToProgUnits(float(self.lambda0),"nm")
            self.use_wavelength = True

            self.Save(use_data_dir=True)
            self.created_rc = True

    def Save(self,use_data_dir=False):
        """Save the preferences"""
        config = ConfigParser.ConfigParser()
        config.add_section("prefs")
#        if not use_data_dir:
#            self.data_dir = self.fcbDataDirectory.get_filename()
        config.set("prefs","data_dir",self.data_dir)

        config.set("prefs","maximize",int(self.maximize))
        config.set("prefs","2n",int(self.power))
        config.set("prefs","local_error",self.local_error)
        config.set("prefs","update_interval",self.update_interval)
        if self.use_wavelength:
            config.set("prefs","graphics","wavelength")
            config.set("prefs","lambda0",self.lambda0)
            config.set("prefs","lambda1",self.lambda1)
        else:
            config.set("prefs","graphics","freq")

        fp = open(self.rc_filename,"w")
        config.write(fp)
        fp.close()

class PreferencesWindow:
    def __init__(self,prefs_):
        global prefs
        prefs = prefs_
        Widgets(os.path.join("ui","preferences.ui")).connect(self)
        self.window = self.winPreferences

        if prefs.created_rc:
            Message(self.window,"Default preferences file created at "+prefs.rc_filename+"\nThe data directory is at "+prefs.data_dir)
        self.UpdateWindow()
        self.fcbDataDirectory.set_current_folder(prefs.data_dir)
#        self.UpdatePrefs()

        self.window.connect("delete_event",self.delete_event)

    def delete_event(self,widget,event,data=None):
        """Do not destroy this window, just hide it"""
        self.window.set_property("visible",False)
        return True

    def on_cmdPrefOk_clicked(self,*args):
        if not self.ValidatePrefs():
            return
        self.UpdatePrefs()
        prefs.Save()
        self.window.set_property("visible",False)

    def on_cmdPrefCancel_clicked(self,*args):
        self.window.set_property("visible",False)

    def LoadPrefs(self,*args):
        """Load the preferences into the controls"""
        self.txt2n.set_text(str(prefs.power))
        self.txtLocalError.set_text(str(prefs.local_error))
        self.txtUpdateInterval.set_text(str(prefs.update_interval))
        self.radioWavelength.set_active(prefs.use_wavelength)

        if prefs.use_wavelength:
            self.txtLambda0.set_text(str(prefs.lambda0))
            self.txtLambda1.set_text(str(prefs.lambda1))

        self.chkMaximize.set_active(prefs.maximize)

    def UpdatePrefs(self,*args):
        """Update the preferences stored on the prefs module"""
        prefs.power = int(self.txt2n.get_text())
        prefs.n = 2**prefs.power
        prefs.local_error = float(self.txtLocalError.get_text())
        prefs.update_interval = float(self.txtUpdateInterval.get_text())
        prefs.use_wavelength = self.radioWavelength.get_active()
        if prefs.use_wavelength:
            prefs.freq_min = units.c/units.ToProgUnits(float(self.txtLambda1.get_text()),"nm")
            prefs.freq_max = units.c/units.ToProgUnits(float(self.txtLambda0.get_text()),"nm")
            prefs.lambda0 = self.txtLambda0.get_text()
            prefs.lambda1 = self.txtLambda1.get_text()

        if self.chkMaximize.get_active():
            prefs.maximize = True
        else:
            prefs.maximize = False

    def ValidatePrefs(self):
        """Ensure the input fields have valid values"""

        if not isFloat(self.txtLocalError.get_text()):
            Error(self.window,"The local error is not valid.")
            return False

        if self.radioWavelength.get_active():
            if not isFloat(self.txtLambda0.get_text()):
                Error(self.window,"The wavelength range start is not valid.")
                return False

            if not isFloat(self.txtLambda1.get_text()):
                Error(self.window,"The wavelength range end is not valid.")
                return False

            if float(self.txtLambda0.get_text())<=0 or float(self.txtLambda1.get_text())<=0:
                Error(self.window,"The wavelengths must be positive.")
                return False

            if float(self.txtLambda0.get_text())>float(self.txtLambda1.get_text()):
                Error(self.window,"The first wavelength must be smaller than the second.")
                return False



        if not isFloat(self.txtUpdateInterval.get_text()):
            Error(self.window,"The update interval is not valid.")
            return False

        if not isInt(self.txt2n.get_text()):
            Error(self.window,"The number of points is not valid.")
            return False

        power = int(self.txt2n.get_text())
        if power<=0:
            Error(self.window,"The number of points is not valid.")
            return False
        elif power>30:
            Error(self.window,"The number of points is excessive, please use a smaller value.")
            return False


        return True

    def show(self):
        self.fcbDataDirectory.set_current_folder(prefs.data_dir)
        self.LoadPrefs()

        self.window.set_property("visible",True)

    def UpdateWindow(self):
        """Enables or disables the widgets according to the selected options"""
        if self.radioFreq.get_active():
            self.hboxWavelengthLimit.set_sensitive(False)
        else:
            self.hboxWavelengthLimit.set_sensitive(True)

    def on_radioFreq_toggled(self,*args):
        self.UpdateWindow()

units = units_module.Units()
