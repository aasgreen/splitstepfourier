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

#import sys


class UnknownUnit(BaseException):
    pass

class UnknownUnitType(BaseException):
    pass

class BasicUnit():
    def __init__(self,type,unit,mult,exp):
        self.type = type
        self.unit= unit
        self.mult = mult
        self.exp = exp

class EquivUnit():
    def __init__(self,unit,unit_equiv,mult,exp):
        self.unit = unit
        self.unit_equiv = unit_equiv
        self.mult = mult
        self.exp = exp

class UnitListItem():
    def __init__(self,unit,exp):
        self.unit = unit
        self.exp = exp
    def __repr__(self):
        return "Unit: %s Exp: %d\n" % (self.unit,self.exp)

#-----------------------------------------------------------------------
# Units class
#-----------------------------------------------------------------------
class Units:
    TYPE_LENGTH,TYPE_TIME,TYPE_CHARGE,TYPE_MASS=range(4)
    desc = ["nm","um","mm","m","km","as","fs","ps","ns","us","ms","s","q","coul","me","g","kg"]
    type = [TYPE_LENGTH,TYPE_LENGTH,TYPE_LENGTH,TYPE_LENGTH,TYPE_LENGTH,TYPE_TIME,TYPE_TIME,TYPE_TIME,TYPE_TIME,TYPE_TIME,TYPE_TIME,TYPE_TIME,TYPE_CHARGE,TYPE_CHARGE,TYPE_MASS,TYPE_MASS,TYPE_MASS]
    mantissa = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.602176462,1.0,9.10938188,1.0,1.0]
    exponent = [-9,-6,-3,0,3,-18,-15,-12,-9,-6,-3,0,-19,0,-31,-3,0]

    basic_units = [
                    # Length
                    BasicUnit('L',"ang",1.0,-10),
                    BasicUnit('L',"nm",1.0,-9),
                    BasicUnit('L',"um",1.0,-6),
                    BasicUnit('L',"mm",1.0,-3),
                    BasicUnit('L',"cm",1.0,-2),
                    BasicUnit('L',"m",1.0,0),
                    BasicUnit('L',"km",1.0,3),
                    # Time
                    BasicUnit('T',"as",1.0,-18),
                    BasicUnit('T',"fs",1.0,-15),
                    BasicUnit('T',"ps",1.0,-12),
                    BasicUnit('T',"ns",1.0,-9),
                    BasicUnit('T',"us",1.0,-6),
                    BasicUnit('T',"ms",1.0,-3),
                    BasicUnit('T',"s",1.0,0),
                    BasicUnit('T',"hr",3.6,3),
                    BasicUnit('T',"day",8.64,4),
                    # Charge
                    BasicUnit('C',"q",1.602176462,-19),
                    BasicUnit('C',"coul",1.0,0),
                    # Mass
                    BasicUnit('M',"me",9.10938188,-31),
                    BasicUnit('M',"g",1.0,-3),
                    BasicUnit('M',"kg",1.0,0)
                  ]

    equiv_units = [
			        EquivUnit("J","kg m^2 s^-2",1.0,0),
			        EquivUnit("uJ","kg m^2 s^-2",1.0,-6),
        			EquivUnit("nJ","kg m^2 s^-2",1.0,-9),
        			EquivUnit("W","kg m^2 s^-3",1.0,0),
        			EquivUnit("MW","kg m^2 s^-3",1.0,6),
        			EquivUnit("GW","kg m^2 s^-3",1.0,9),
        			EquivUnit("eV","kg m^2 s^-2",1.60217733,-19),
        			EquivUnit("V","kg m^2 s^-2 coul^-1",1.0,0)
                  ]

    _instance = None
    def __new__(cls, *args, **kwargs):  #Singleton (Class with a single instance)

        if not cls._instance:
            cls._instance = super(Units, cls).__new__(cls, *args, **kwargs)
        return cls._instance

    def __init__(self):
        self.SetProgramUnits("um","fs","me","q")
#        print "1 W",self.ToProgUnits(1.0,"W")
#        self.SetProgramUnits("um","fs","me","q")
#        print "1E-12 s^-1",self.FromProgUnits(1.0E-12,"s^-1")
#        print "1 W",self.ToProgUnits(1.0,"W")
#        print "c",c,"SI:",self.FromProgUnits(c,"m s^-1")
#        sys.exit(1)

    def DefineConstants(self):
        self.c = self.ToProgUnits(299792458.0,"m s^-1")
        self.epson0 = self.ToProgUnits(8.85418781762E-12,"coul^2 s^2 m^-3 kg^-1")
        self.me = self.ToProgUnits(9.10938188E-31,"kg")
        self.q = self.ToProgUnits(1.602176462E-19,"coul")
        self.hbar = self.ToProgUnits(1.05457266913E-34,"J s")

    def SetProgramUnits(self,length_unit,time_unit,mass_unit,charge_unit):
        """Sets the program units. Arguments are strings"""

        self.length_unit = length_unit
        self.time_unit = time_unit
        self.mass_unit = mass_unit
        self.charge_unit = charge_unit
        self.L = self.index(Units.basic_units,self.length_unit)
        self.T = self.index(Units.basic_units,self.time_unit)
        self.M = self.index(Units.basic_units,self.mass_unit)
        self.C = self.index(Units.basic_units,self.charge_unit)

        self.DefineConstants()

    def index(self,list,unit):
        i = 0
        for e in list:
            if e.unit == unit:
                return i
            i += 1

        return -1

    def ToProgUnits(self,value,from_unit):
        return self.Convert(value,from_unit,True)

    def FromProgUnits(self,value,from_unit):
        return self.Convert(value,from_unit,False)

    # -----------------------------------------------------------------------
    # Needs expansion?
    # -----------------------------------------------------------------------
    def NeedsExpansion(self,unit_list):
        for e in unit_list:
            if self.index(Units.equiv_units,e.unit)!=-1:
                return True
        return False

    # -----------------------------------------------------------------------
    # Replace the element pos of unit_list by equiv_list
    # -----------------------------------------------------------------------
    def Replace(self,unit_list,equiv_list,pos):
        return unit_list[:pos]+equiv_list+unit_list[pos+1:]

    # -----------------------------------------------------------------------
    # Expand
    # -----------------------------------------------------------------------
    def Expand(self,unit_list):
        if not self.NeedsExpansion(unit_list):
            return (unit_list,1.0,0)

        mult = 1.0
        exp = 0

        n = len(unit_list)

        for i in range(n):
            eq_ind = self.index(Units.equiv_units,unit_list[i].unit)
            if eq_ind != -1:
                equiv_list = self.Tokenize(Units.equiv_units[eq_ind].unit_equiv)
                m = len(equiv_list)
                for j in range(m):
                    equiv_list[j].exp *= unit_list[i].exp

                unit_list = self.Replace(unit_list,equiv_list,i)
                mult *= Units.equiv_units[eq_ind].mult**unit_list[i].exp
                exp += Units.equiv_units[eq_ind].exp*unit_list[i].exp

        return (unit_list,mult,exp)

    # -----------------------------------------------------------------------
    # Tokenize
    # -----------------------------------------------------------------------
    def Tokenize(self,value_unit):
        unit_list = []

        u = value_unit.split(" ")

        exp_coeff = 0
        for e in u:
            exp_pos = e.find("^")
            if exp_pos != -1:
                l = e.split("^")
                unit = l[0]
                exp = int(l[1]) #FIXME: Add support to fractional exponents?
            else:
                unit = e
                exp = 1

            unit_list.append(UnitListItem(unit,exp))

        return unit_list


    # -----------------------------------------------------------------------
    # Returns the program unit index for this type
    # -----------------------------------------------------------------------
    def ProgUnitInd(self,type):
        if type == "L":
            return self.L
        elif type == "T":
            return self.T
        elif type == "M":
            return self.M
        elif type == "C":
            return self.C
        else:
            print "ProgUnitInd: Unknown type",type
            raise UnknownUnitType

    def Convert(self,value,value_unit,to_prog_units):
        """Convert value to program units.
        value may be a list or an array.
        from_unit: The various entries are separated by spaces, "^" means exponentiation
        Example: ToProgUnits(1.0,"m s^-1")
        """

        unit_list = self.Tokenize(value_unit)
        unit_list,expand_mult,expand_exp = self.Expand(unit_list)
        
        mult = 1.0
        exp = 0

        n = len(unit_list)

        for i in range(n):
            l_ind = self.index(Units.basic_units,unit_list[i].unit)
            if l_ind == -1:
                print "Unknown unit",unit_list[i].unit
                raise UnknownUnit
            l_exp = unit_list[i].exp
            p_ind = self.ProgUnitInd(Units.basic_units[l_ind].type)
            if to_prog_units:
                exp += (Units.basic_units[l_ind].exp-Units.basic_units[p_ind].exp)*l_exp
                mult *= (Units.basic_units[l_ind].mult/Units.basic_units[p_ind].mult)**l_exp
            else:
                exp += (Units.basic_units[p_ind].exp-Units.basic_units[l_ind].exp)*l_exp
                mult *= pow((Units.basic_units[p_ind].mult/Units.basic_units[l_ind].mult),l_exp)

        if to_prog_units:
            return value*mult*(10.0**(exp+expand_exp))*expand_mult
        else:
            return value*mult*(10.0**(exp-expand_exp))/expand_mult

