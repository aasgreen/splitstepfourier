# -*- coding: latin1 -*-
#  Copyright (C) 2009 João Luís Silva <jsilva@fc.up.pt>
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

import os
import sys

if hasattr(sys,"frozen"):
    basedir = os.path.split(os.path.dirname(__file__))[0]
else:
    basedir = os.path.dirname(__file__)


class Widgets:
    def __init__(self,file):
        full_path_to_file = os.path.join(basedir,file)
        self.widgets = gtk.Builder()
        self.widgets.add_from_file(full_path_to_file)

    def connect(self,obj,signals_obj=None):
        self.connect_items(obj)
        if signals_obj:
            self.connect_signals(signals_obj)
        else:
            self.connect_signals(obj)

    def connect_items(self,obj):
        """Creates new variables bound to the corresponding widgets"""
        l = self.widgets.get_objects()
        for e in l:
            # Using gtk.Buildable due to a bug in GTK, see https://bugs.launchpad.net/ubuntu/+source/pygtk/+bug/507739
            name = gtk.Buildable.get_name(e)    
            setattr(obj, name, self[name])

    def connect_signals(self,obj):
        """Connects all obj methods starting with on_ and after_"""
        self.widgets.connect_signals(obj)

    def __getitem__(self,key):
        return self.widgets.get_object(key)
