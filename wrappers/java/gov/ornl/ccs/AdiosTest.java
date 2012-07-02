package gov.ornl.ccs;

import gov.ornl.ccs.Adios;
import gov.ornl.ccs.AdiosFile;
import gov.ornl.ccs.AdiosGroup;

public class AdiosTest
{
    // The main program
    public static void main(String[] args)
    {
        System.out.println("Hello");

        //Adios adios = new Adios();
        //adios.init();

        AdiosFile file = new AdiosFile();
        file.open("adios_global.bp");

        AdiosGroup group = new AdiosGroup(file);
        group.open("temperature");

        AdiosVarinfo var = new AdiosVarinfo(group);
        var.inq("temperature");

        long[] start = {0, 0};
        long[] count = {1, 10};
        double[] output = var.read(start, count);
        //double[] output = var.read(start, count);

        System.out.println("output.length = " + output.length);
        for (int i = 0; i < output.length; i++) {
            System.out.println("output[" + i + "] = " + output[i]);
        }

        var.close();
        group.close();
        file.close();

        //AdiosGroup group = new AdiosGroup();
        //group.open("adios_global.bp");
        //group.open("temperature");

        
        System.out.println();
        System.out.println(file);
        System.out.println(group);
        System.out.println(var);
    }
}
