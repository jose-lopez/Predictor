 /*
 GenInformation.java
 Esta clase es la que construtuye el archivo de lecturas utilizando la clase Middle, esta a su vez es la que conecta con prolog--------

 Copyright (C) 2015 
 Morales Yonathan (yonathan.morales@unet.edu.ve), Jose Lopez (jlopez@unet.edu.ve).

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>

 */
// Esta clase es la que construtuye el archivo de lecturas utilizando la clase Middle, esta a su vez es la que conecta con prolog--------
package util;

import Libreria.Utilities;
import gene.feature.Exon;
import gene.feature.Gene;
import gene.feature.Information;
import gene.feature.UTR3p;
import gene.feature.UTR5p;
import gene.information.Analizer;
import gene.information.GeneConstructor;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import pipeline.BioPattern;
import pipeline.Region;
import util.MiddleWare;

public class GenInformation {
    //-- estas 4 variables boolean se utilizan para excluir lo que se va a mostrar en el archivo final, que son los atg, las paradas, los intrones y exones----

    private boolean STOP;
    private boolean INTRONES;
    private boolean EXONES;
    private String hebra;
    private String lect, KT;
    private int[] globIDs = {0, 0};
    private boolean ATG;
    private boolean GORF;
    private boolean GT;
    private boolean GPR;
    private boolean complementarios;

    //--- este constructor recibe la hebra, el cual pone todo en true para que haga todo y lect completo para que haga todas las lecturas---
    public GenInformation(String Hebra) {
        ATG = true;
        STOP = true;
        INTRONES = true;
        EXONES = true;
        hebra = Hebra;
        lect = "completo";
    }
//Este constructor recibe la hebra, recibe las lecturas y los boolean de los atg paradas ect , para especificar que es lo que se quiere hacer-----------

    public GenInformation(String Hebra, String Lecturas, boolean atg, boolean stop, boolean intrones, boolean exones, boolean gorf, boolean gt, boolean gpr, boolean complementarios) {
        ATG = atg;
        STOP = stop;
        INTRONES = intrones;
        EXONES = exones;
        hebra = Hebra;
        lect = Lecturas;
        this.setGORF(gorf);
        this.setGT(gt);
        this.setGPR(gpr);
        this.setComplementarios(complementarios);
    }
    //--genera recibe el archivo temporal creado en la funcion inicio ----

    public boolean isComplementarios() {
        return complementarios;
    }

    public void setComplementarios(boolean complementarios) {
        this.complementarios = complementarios;
    }

    public boolean isGORF() {
        return GORF;
    }

    public boolean isGT() {
        return GT;
    }

    public boolean isGPR() {
        return GPR;
    }

    public void setGORF(boolean GORF) {
        this.GORF = GORF;
    }

    public void setGT(boolean GT) {
        this.GT = GT;
    }

    public void setGPR(boolean GPR) {
        this.GPR = GPR;
    }

    // Probando probando
    public int[] getGlobIDs() {
        return globIDs;
    }

    public boolean isATG() {
        return ATG;
    }

    public boolean isSTOP() {
        return STOP;
    }

    public boolean isINTRONES() {
        return INTRONES;
    }

    public boolean isEXONES() {
        return EXONES;
    }

    public String getHebra() {
        return hebra;
    }

    public String getLect() {
        return lect;
    }

    public String getKT() {
        return KT;
    }

    public void guardar_2(String palabra, File salida) throws IOException {
        try (BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(salida, false), "UTF8"))) {
            out.write(palabra);
            out.write("\n");
        }
    }
//-- esta funcion es llamada desde principal la cual envia los url de los archivos 
// genesEnProceso, idsGenesEnProceso, salidaPredGFF3, salidaEnsEPDGFF3

    public void inicio(String genesEnProceso, String idsGenesEnProceso, String salidaPredGFF3, String salidaEnsEPDGFF3, String gffVersion, String regionID, String regionInicio, String regionFin, boolean ilpinr, boolean consensos, boolean reporteAbs, int numObjs, int numIter, boolean ilpClasificador, String red) throws IOException, Exception {
        //--se usa el archivo de genesEnProceso y se le pasa a la clase utiliti de la libreria-----------       
        Utilities auxRegistrosEnsembl = new Utilities(genesEnProceso, hebra);
        //-- se desglosa el archivo de genesEnProceso y se escribe un archivo temporal,se hace un archivo temporal por cada gen q se encuentre en la genesEnProceso-------       
        String escritura = "";
        File au = new File("tem");
        File gff3Predictor = new File(salidaPredGFF3);
        File gff3EnsEPD = new File(salidaEnsEPDGFF3);
        gff3Predictor.delete();
        gff3EnsEPD.delete();

        // Se imprime encabezado archivo de lecturas en el archivo salida.gff3Predictor
        String versionGFF = "##gff-version\t" + gffVersion + "\n";
        String encabRegion = "##sequence-region\t" + regionID + "\t" + regionInicio + "\t" + regionFin + "\n";

        auxRegistrosEnsembl.guardar(versionGFF + encabRegion, gff3Predictor);
        auxRegistrosEnsembl.guardar(versionGFF + encabRegion, gff3EnsEPD);


        int genIDsize = auxRegistrosEnsembl.get_GenID().size();
        for (int i = 0; i < genIDsize; i++) {
            au.delete();

            escritura = "(["
                    + auxRegistrosEnsembl.get_GenID().get(i) //  Gen ID del Archivo Original
                    + "],[" + auxRegistrosEnsembl.get_InicioLocal().get(i) //  Inicio Local del Archivo Original  
                    + "],[" + auxRegistrosEnsembl.get_FinLocal().get(i) //  Fin Local del Archivo Original
                    + "]," + auxRegistrosEnsembl.get_Secuencia().get(i) //  Secuencia del Gen
                    + ",[" + auxRegistrosEnsembl.get_Cromosoma().get(i) //  Cromosoma
                    + "]," + auxRegistrosEnsembl.get_InicioAbsoluto().get(i) //  Inicio Absoluto del Archivo Original
                    + "," + auxRegistrosEnsembl.get_FinAbsoluto().get(i) //  Fin Absoluto del Archivo Original
                    + ")"
                    + "([0,0],[0,0],[0])";

            this.guardar_2(escritura, au);
            //---llama a la funcion genera que es la que genera las lecturas en el archivo de salida-------------          
            this.generaLects("tem", gff3EnsEPD, gff3Predictor, ilpinr, consensos, reporteAbs, numObjs, numIter, ilpClasificador, red);
            au.delete();
        }

    }

    public void setGlobIDs(int[] globIDs) {
        this.globIDs = globIDs;
    }

    public void setATG(boolean ATG) {
        this.ATG = ATG;
    }

    public void setSTOP(boolean STOP) {
        this.STOP = STOP;
    }

    public void setINTRONES(boolean INTRONES) {
        this.INTRONES = INTRONES;
    }

    public void setEXONES(boolean EXONES) {
        this.EXONES = EXONES;
    }

    public void setHebra(String hebra) {
        this.hebra = hebra;
    }

    public void setLect(String lect) {
        this.lect = lect;
    }

    public void setKT(String KT) {
        this.KT = KT;
    }

    public void generaLects(String entrada, File gff3EnsemblEPDExt, File gff3Predictor, boolean iLPinr, boolean consensos, boolean reporteAbs, int numObjs, int numIter, boolean ilpClasificador, String red) throws IOException, Exception {

        /* Coordenadas VEGA y Ensembl SST
         List<Integer> atg = new ArrayList<>(Arrays.asList(108));
         List<Integer> gt = new ArrayList<>(Arrays.asList(246));
         List<Integer> ag = new ArrayList<>(Arrays.asList(1121));
         List<Integer> stops = new ArrayList<>(Arrays.asList(1334));
         //*/
        
         /* Coordenadas VEGA y Ensembl SST + 2000 up and down
         List<Integer> atg = new ArrayList<>(Arrays.asList(2108));
         List<Integer> gt = new ArrayList<>(Arrays.asList(2246));
         List<Integer> ag = new ArrayList<>(Arrays.asList(3121));
         List<Integer> stops = new ArrayList<>(Arrays.asList(3334));
         //*/

         /* Coordenadas vector juguete
         List<Integer> atg = new ArrayList<>(Arrays.asList(5));
         List<Integer> gt = new ArrayList<>(Arrays.asList(23));
         List<Integer> ag = new ArrayList<>(Arrays.asList(38));
         List<Integer> stops = new ArrayList<>(Arrays.asList(56));
         //*/

        /*
         List<Integer> atg = new ArrayList<>(Arrays.asList(1721,2101,2749));
         List<Integer> gt = new ArrayList<>(Arrays.asList(1818,2239,2541));
         List<Integer> ag = new ArrayList<>(Arrays.asList(2655,3114,3141));
         List<Integer> stops = new ArrayList<>(Arrays.asList(1536,3327,3569));
         //*/

        /*
         List<Integer> atg = new ArrayList<>(Arrays.asList(1721, 2101));
         List<Integer> gt = new ArrayList<>(Arrays.asList(1818, 2239));
         List<Integer> ag = new ArrayList<>(Arrays.asList(2655, 3114));
         List<Integer> stops = new ArrayList<>(Arrays.asList(1536, 3327));//*/

        /*
         List<Integer> atg = new ArrayList<>(Arrays.asList(2101));
         List<Integer> gt = new ArrayList<>(Arrays.asList(1847, 2239, 2459));
         List<Integer> ag = new ArrayList<>(Arrays.asList(2059, 2341, 2655, 3114));
         List<Integer> stops = new ArrayList<>(Arrays.asList(3327));//*/


        //---------------------Nombre de los archivos------------------------------------------------
        String url_archivo_entrada = entrada; // url del Archivo de entrada de la base de datos de internet se puede usar args[0]
        String url_archivo_middle = "gen_prueba.pl";     // url gen problema que se usa para la clase Middle
        String rutaGen = "salidas/gen.txt";     // url gen problema que se usa para la clase Middle
        String rutaGenClasificador = "salidas/genCla.txt";     // url gen problema que se usa para la clase Middle
        //------------------Archivos-------------------------------------------------------------------

        File genes = new File(url_archivo_middle);// archivo .pl de entrada a la clase middle
        File gen = new File(rutaGen);// archivo .pl de entrada a la clase middle
        File genDos = new File(rutaGenClasificador);// archivo .pl de entrada a la clase middle

        Utilities metaDataGen = new Utilities(url_archivo_entrada, hebra);
        //--------------------Codigo para generar el archivo de entrada al Middle-----------------------
        String secuenciaProceso = metaDataGen.get_Secuencia().get(0);
        String aux = "gen(" + secuenciaProceso + ")."; //se genera gen_prueba
        //String aux2 = secuenciaProceso; //se genera gen para clasificador

        String secuencia = metaDataGen.get_Secuencia().toString(), secuenciaDos;
        secuencia = secuencia.replaceAll("[\\[\\]]", "");
        secuenciaDos = secuencia.replaceAll("[\\[\\]]", "");
        secuencia = secuencia.replaceAll(",", "");
        //.replaceAll("[", "").replaceAll("]", "");

        genes.delete();
        gen.delete();
        genDos.delete();

        this.guardar_2(aux, genes);
        this.guardar_2(secuencia, gen);
        this.guardar_2(secuenciaDos, genDos);

        List<String> data;
        data = readGene(rutaGen);

        MiddleWare middle = new MiddleWare(); // este es el objero de clase midddleware que conecta con el predictor

        middle.init("p_genes.pl"); //Se inicializa el objeto predictor con el codigo del predictor.

        Analizer analizer = new Analizer();
        //analizer.readFromLists(atg, gt, ag, stops, data); // Descomentar para trabajar con listas

        analizer.readFromMiddleWare(middle, ilpClasificador, data, rutaGenClasificador, secuencia); // Se instancia el objeto constructor de lecturas
        // empleando las predicciones desde Prolog que estan disponibles en el objeto middle. 
        // Al salir de este metodo las listas de predicciones hechas desde prolog estan disponibles
        // en las listas atg, gt, ... del objeto constructor presente en el objeto analizer.
        // Descomentar para trabajar con prolog.*/
        GeneConstructor constructor = analizer.getConstructor(); // Se puede acceder a las predicciones.

        analizer.constructLectures(false); //Se construyen las lecturas (true metodo recursivo, false interativo).

        System.out.println("" + analizer.toString()); // Se imprimen en consola todos los intrones y exones hallados.

        if (!analizer.getLectures().isEmpty()) {

            System.out.println("Numero de ORFs: " + analizer.getLectures().size());

            if (this.isGT()) {
                analizer.constructRegionUTR5p(metaDataGen, this, iLPinr, consensos); // pasar iLPinr = false implica que se proponen TSSs por consensos.
            }
            //analizer.constructRegionUTR3p(metaDataGen, this, false);
            
            analizer.lectsEnsemblEpd(metaDataGen, gff3EnsemblEPDExt, reporteAbs, this);

            analizer.lectsToGFF3(metaDataGen, gff3Predictor, reporteAbs, this, numObjs, numIter, red, this.isComplementarios());

            metaDataGen.guardar("\n", gff3Predictor);

        } else {
            System.out.println("No hay lecturas para " + metaDataGen.get_GenID().get(0));

        }
    }
//-- genera posibles lecturas en la region UTR5'. Necesaria para reportar posibles exones fuera del ORF.

    public Analizer generaLectsUTR5p(List<Information> regionUTR5) throws IOException, Exception {

        String url_archivo_middle = "gen_prueba.pl";     // url gen problema que se usa para la clase Middle
        //String url_archivo_inter = URL_Intermedio; // url Salida intermedia usada para generar el GTF
        //String url_archivo = "ghy.txt";
        //------------------Archivos-------------------------------------------------------------------
        //File temp = new File("temporal.txt");// archivo temporal de las lecturas
        File utr5p = new File(url_archivo_middle);// archivo .pl de entrada a la clase middle
        //File inter = new File(url_archivo_inter);  //    archivo intermedio con las lecturas que da el middleware
        //File au = new File(url_archivo);
        // se le envia a la clase utilities de libreria GTF el archivo de entrada para extraer todos los datos de los utr5p
        //--------------------Codigo para generar el archivo de entrada al Middle-----------------------
        ArrayList<String> utr5 = new ArrayList<>();
        for (Information i : regionUTR5) {
            utr5.add(i.toString());
            /*if (i.position < regionUTR5.size()) {
             utr5.add(",");
             }*/
        }

        String aux = utr5.toString().replaceAll(" ", "");
        aux = "gen(" + aux + ")."; //se genera gen_prueba

        utr5p.delete();
        this.guardar_2(aux, utr5p);
        //metaDataGen.set_SecLect(lect);//se indica a Utilities cuales lecturas del gen se reportaran en el GTF.

        //KT = " ";
        // comentar para desactivar el middleware y trabajar solo con las listas , 
        // ya que este programa se puede trabajar con listas o con  Middle por lo tanto con prolog.
        //*

        MiddleWare middle = new MiddleWare(); // este es el objero de clase midddleware que conecta con el predictor
        middle.init("p_genes.pl"); //Se inicializa el objeto predictor con el codigo del predcitor.
        //middle.consultEverything(); //hacen las consultas y se cargan los vectores Patgs, Pgt, Pag, Ppar.
        //*/

        Analizer analizer = new Analizer();

        analizer.readFromMiddleWare(middle, true, utr5, aux, aux); // Se instancia el objeto constructor de lecturas
        // empleando las predicciones desde Prolog que estan disponibles en el objeto middle. 
        // Al salir de este metodo las listas de predicciones hechas desde prolog estan disponibles
        // en las listas atg, gt, ... del objeto constructor presente en el objeto analizer.
        // Descomentar para trabajar con prolog.*/


        GeneConstructor constructor = analizer.getConstructor();
        //constructor.setAtg(new ArrayList<Integer>());
        //constructor.setStops(new ArrayList<Integer>());


        //analizer.readFromLists(atg, gt, ag, stops, data); // Descomentar para trabajar con listas
        analizer.constructLectures(false); //Se construyen las lecturas (true metodo recursivo, false interactivo).
        //System.out.println("" + analizer.toString()); // Se imprimen en consola todos los intrones y exones hallados.

        //analizer.filtrar();
        //analizer.lectsToString();


        return analizer;

    }

    public static List readGene(String rutaGen) throws FileNotFoundException, IOException {
        String path = System.getProperty("user.dir");
        path += "/" + rutaGen;
        //Scanner s = new Scanner(new File(path));
        BufferedReader br = new BufferedReader(new FileReader(path));
        ArrayList<String> list = new ArrayList<String>();
        int r;
        while ((r = br.read()) != -1) {
            char ch = (char) r;
            list.add(Character.toString(ch).toLowerCase());
        }
        list.remove(list.size() - 1);
        br.close();

        return list;
    }
}
