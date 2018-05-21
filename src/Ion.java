//package simSpec;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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