//package simSpec;

/*
 *  Copyright, 2018, Muaaz Awan and Fahad Saeed
    This file is part of MaSS-Simulator

    MaSS-Simulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MaSS-Simulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MaSS-Simulator.  If not, see <https://www.gnu.org/licenses/>.
 *
 *
 */
//package specsimulate;

/**
 *
 * @author Gul
 */
public class Ion {
    public double m_z;
    public double intensity;
    public Ion(){
        m_z = 0;
        intensity = 0;
    }
    
    public Ion(double m_zN, double intensityN){
        m_z = m_zN;
        intensity = intensityN;
    }
    
    
    int compareTo(Ion myIon){
        if(this.m_z > myIon.m_z)
            return 1;
        else if(this.m_z < myIon.m_z)
            return -1;
        else
            return 0;
    }
}
