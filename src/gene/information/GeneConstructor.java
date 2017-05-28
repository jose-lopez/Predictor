 /*
 GeneConstructor.java
 Clase para manegar la construccion de genes(lecturas) posibles usando las
 listas de inicios(atg), paradas(stops), y intrones (gt,ag) como indices
 referentes a sus valores en la lista de informacion(geneData)

 Copyright (C) 2016. 
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
package gene.information;

import gene.feature.Information;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import util.MiddleWare;

/**
 * Clase para manegar la construccion de genes(lecturas) posibles usando las
 * listas de inicios(atg), paradas(stops), y intrones (gt,ag) como indices
 * referentes a sus valores en la lista de informacion(geneData)
 */
public class GeneConstructor {
    //---------------------------Private Attributes-----------------------------
    // <editor-fold desc="Private Attributes">
//---------- listas de inicio, como indices referentes a la informacion que se encuentra en gendata---------------------------------

    private List<Integer> atg;
    private List<Integer> gt;
    private List<Integer> ag;
    private List<Integer> stops; //Pueden ser: taa, tag o tga. Se guardan sus posiciones.
    //---------------------------------------
    private List<Information> geneData;
    //---------------------------------------
    private boolean withoutStarts = false;
    private boolean withoutStops = false;
    //  </editor-fold>    
    //---------------------------Setters---------------------------------------- 
    // <editor-fold defaultstate="collapsed" desc="Setters">

    public void setAtg(List<Integer> atg) {
        this.atg = atg;
    }

    public void setGt(List<Integer> gt) {
        this.gt = gt;
    }

    public void setAg(List<Integer> ag) {
        this.ag = ag;
    }

    public void setStops(List<Integer> stops) {
        this.stops = stops;
    }

    public void setGeneData(List<Information> geneData) {
        this.geneData = geneData;
    }

    public void setWithoutStarts(boolean withoutStarts) {
        this.withoutStarts = withoutStarts;
    }

    public void setWithoutStops(boolean withoutStops) {
        this.withoutStops = withoutStops;
    }

    //  </editor-fold>
    //---------------------------Constructors-----------------------------------
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor usado para crear el objeto con las listas ya llenas de
     * acuerdo a los valores que llegan por parametros
     */
    public GeneConstructor(List<Integer> atg, List<Integer> gt, List<Integer> ag, List<Integer> stops, List<Information> geneData) {
        this.atg = atg;
        this.gt = gt;
        this.ag = ag;
        this.stops = stops;
        this.geneData = geneData;

        isCompatibleGene();
    }

    //---------------------------------------
    /**
     * Constructor usado para llenar las listas de uso interno a traves de un
     * objeto de la clase Middle ya instanciado y probado, de donde solo se
     * obtendran los valores que este devuelve como entrada para generar las
     * listas internas
     */
    public GeneConstructor(MiddleWare middleWare, boolean ilpClasificador, List<String> data , String rutaSecuencia,String secuencia) throws Exception {

        if (ilpClasificador) {
            this.initLists(middleWare.getAtgPositions(),
                    middleWare.getGtPositions(),
                    middleWare.getAgPositions(),
                    middleWare.getParadasPositions(),
                    middleWare.getGenData());
        }else{
            
            
            List<Integer> gts =  middleWare.getGtPositionsClasificador(0, 1, rutaSecuencia);
            List<Integer> ags =  middleWare.getAGPositionsClasificador(1, 1, rutaSecuencia);
            List<Integer> atgs = middleWare.getPositionsPatron(secuencia, true);
            List<Integer> stopss = middleWare.getPositionsPatron(secuencia, false);
            this.initLists(atgs, gts, ags, stopss, middleWare.getGenData());
            
        }

        this.isCompatibleGene();
    }

    //---------------------------------------
    /**
     * Constructor usado para llenar las listas de uso interno a traves de un
     * objeto de la clase Middle ya instanciado y probado, de donde solo se
     * obtendran los valores que este devuelve como entrada para generar las
     * listas internas
     */
    /*public GeneConstructor(Middle middleWare) throws Exception {
     this.initLists(new ArrayList<Integer>(),
     middleWare.getGtPositions(),
     middleWare.getAgPositions(),
     new ArrayList<Integer>(),
     middleWare.getGenData());

     this.isCompatibleGene();
     }*/
    //---------------------------------------
    /**
     * Constructor solo debe ser usado si el objeto que se esta creando esta en
     * otra clase del mismo paquete pues, el objetivo seria poder instanciar
     * luego todas sus listas a traves del metodo initLists
     */
    public GeneConstructor() {
    }

    //  </editor-fold>
    //---------------------------Getters---------------------------------------- 
    // <editor-fold defaultstate="collapsed" desc="Getters">
    public List<Integer> getAtg() {
        return atg;
    }

    //---------------------------------------
    public List<Integer> getGt() {
        return gt;
    }

    //---------------------------------------
    public List<Integer> getAg() {
        return ag;
    }

    //---------------------------------------
    public List<Integer> getStops() {
        return stops;
    }

    //---------------------------------------
    public int lastData() {
        return geneData.size() - 1;
    }

    //---------------------------------------
    public Information getData(int i) {
        return this.geneData.get(i);
    }

    //---------------------------------------
    public List<Information> getGeneData() {
        return geneData;
    }

    //---------------------------------------
    public List<Information> getInnerInfo(int from, int to) {
        return geneData.subList(from, to);
    }

    //---------------------------------------
    public boolean isWithoutStops() {

        return withoutStops;
    }

    //---------------------------------------
    public boolean isWithoutStarts() {

        return withoutStarts;
    }

    //  </editor-fold>
    //---------------------------Public Methods--------------------------------- 
    // <editor-fold defaultstate="collapsed" desc="Public Methods">
    public static void main(String[] args) throws Exception {

        GeneConstructor constructor = new GeneConstructor();
        int limiteP = 250;
        ArrayList<Integer> atgs = new ArrayList<>();
        for (int i = 0; i < 4; i++) {
            Integer atg = i * 100;
            atgs.add(atg);
        }
        constructor.setAtg(atgs);
        constructor.resetATG(limiteP, "Gen000");

    }

    /**
     * Verificaciones de casos especificos para la compatibilidad de un gen
     */
    public boolean isCompatibleGene() {

        if ((this.atg.isEmpty() || this.stops.isEmpty()) && (this.gt.isEmpty() || this.ag.isEmpty())) {

            System.err.println("ADVERTENCIA: No se pueden armar lecturas"
                    + "\n El predictor no hallo' suficientes coordenadas");
            return false;
        }

        if (!this.gt.isEmpty() && !this.ag.isEmpty()) {
            if (this.atg.isEmpty()) {
                System.err.println("ADVERTENCIA: GEN SIN INICIO, lecturas de marco abierto"
                        + "\n validaciones de ATGs desactivadas");
                this.withoutStarts = true;
                this.atg.add(new Integer(0));
            }
            if (this.stops.isEmpty()) {
                System.err.println("ADVERTENCIA: GEN SIN FIN, lecturas de marco abierto"
                        + "\n validaciones de Paradas desactivadas");
                this.withoutStops = true;
                this.stops.add(new Integer(geneData.size() - 1));
            }
        }

        return true;
    }

    /**
     * Se redefinen los sitios ATG segun las coordenadas del promotor de la
     * secuencia.
     */
    public boolean resetATG(int limiteP, String genID) {

        boolean reset = false;
        int limitATGs;
        if (limiteP > 200) {
            limitATGs = limiteP - 200;
        } else {
            limitATGs = limiteP;
        }

        ArrayList<Integer> auxATG = new ArrayList<>();

        for (Integer a : this.atg) {

            if (!(a < limitATGs)) {
                auxATG.add(a);
            }

        }
        if (auxATG.isEmpty()) {
            reset = false;
            System.out.println("Region promotora solapada con sitios de inicio en: " + genID);
        } else {
            this.atg = auxATG;
            System.out.println("Sitios de inicio reestructurados segun region promotora en: " + genID);
            reset = true;
        }

        return reset;

    }

    //---------------------------------------
    //  </editor-fold>
    //---------------------------Package Methods-------------------------------- 
    // <editor-fold defaultstate="collapsed" desc="Package Methods">
    /**
     * Metodo usado para instanciar las listas a traves de listas de Integers
     * para los indices, y una lista de Strings para la data, solo es posible
     * usarlo si se esta en el mismo paquete que esta clase, por lo que se
     * recomienda usar en lugar de esto los constructores
     */
    void initLists(List<Integer> atg, List<Integer> gt, List<Integer> ag, List<Integer> stops, List<String> data) throws Exception {
        this.atg = atg;
        this.gt = gt;
        this.ag = ag;
        //this.gt = new ArrayList<>();
        //this.ag = new ArrayList<>();
        //this.stops = new ArrayList<>();
        //this.atg = new ArrayList<>();
        this.stops = stops;
        //int temp = 1126;
        //stops.add(temp);
        //this.atg.add(646);
        //Collections.sort(this.stops);
        this.geneData = new ArrayList<>();

        for (int i = 0; i < data.size(); i++) {
            Information info = new Information(new Integer(i), data.get(i));
            this.geneData.add(info);
        }
    }
    //---------------------------------------
    //  </editor-fold>
    //---------------------------Override Methods------------------------------- 
    // <editor-fold defaultstate="collapsed" desc="Override Methods">

    @Override
    public String toString() {
        String out = "";
        String aux = "";

        /*aux = "ATG\t= [";
         for (Integer index : atg) {
         aux += ", " + index;
         }
         aux = aux.replaceFirst(", ", "");
         aux += "]";
        
         out += aux + "\n";
        
         aux = "GT\t= [";
         for (Integer index : gt) {
         aux += ", " + index;
         }
         aux = aux.replaceFirst(", ", "");
         aux += "]";
        
         out += aux + "\n";
        
         aux = "AG\t= [";
         for (Integer index : ag) {
         aux += ", " + index;
         }
         aux = aux.replaceFirst(", ", "");
         aux += "]";
        
         out += aux + "\n";
        
         aux = "STOPS\t= [";
         for (Integer index : stops) {
         aux += ", " + index;
         }*/
        aux = aux.replaceFirst(", ", "");
        aux += "]";

        out += aux + "\n";
        return out;
    }
    //  </editor-fold>
}
