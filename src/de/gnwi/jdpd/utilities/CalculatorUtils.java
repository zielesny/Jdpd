package de.gnwi.jdpd.utilities;

import java.util.concurrent.ExecutorService;

/**
 * Utility class for calculators
 * 
 * @author Achim Zielesny
 */
public class CalculatorUtils {

    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Executor service
     */
    private ExecutorService executorService;
    
    /**
     * Adder groups
     */
    private AdderGroup[] adderGroups;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    public CalculatorUtils() {
        this.executorService = null;
        this.adderGroups = null;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Resets adders
     * NOTE: NO checks are performed.
     */
    public void resetAdders() {
        if (this.adderGroups != null) {
            for (AdderGroup tmpAdderGroup : this.adderGroups) {
                tmpAdderGroup.reset();
            }
        }
    }

    /**
     * Executor service shutdown
     */
    public void shutdownExecutorService() {
        if (this.executorService != null) {
            this.executorService.shutdown();
            this.executorService = null;
        }
    }
    
    /**
     * Return accumulated (total) sum of all potential energy adders
     * NOTE: NO checks are performed.
     * 
     * @return Return accumulated (total) sum of all potential energy adders
     */
    public double getAccumulatedPotentialEnergyAddersSum() {
        double tmpTotalSum = 0.0;
        for (AdderGroup tmpAdderGroup : this.adderGroups) {
            tmpTotalSum += tmpAdderGroup.getPotentialEnergyAdder().getSum();
        }
        return tmpTotalSum;
    }
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal x term adders
     * NOTE: NO checks are performed.
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal x term adders
     */
    public double getAccumulatedPressureXAddersSum() {
        double tmpTotalSum = 0.0;
        for (AdderGroup tmpAdderGroup : this.adderGroups) {
            tmpTotalSum += tmpAdderGroup.getPressureXAdder().getSum();
        }
        return tmpTotalSum;
    }
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal y term adders
     * NOTE: NO checks are performed.
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal y term adders
     */
    public double getAccumulatedPressureYAddersSum() {
        double tmpTotalSum = 0.0;
        for (AdderGroup tmpAdderGroup : this.adderGroups) {
            tmpTotalSum += tmpAdderGroup.getPressureYAdder().getSum();
        }
        return tmpTotalSum;
    }
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal z term adders
     * NOTE: NO checks are performed.
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal z term adders
     */
    public double getAccumulatedPressureZAddersSum() {
        double tmpTotalSum = 0.0;
        for (AdderGroup tmpAdderGroup : this.adderGroups) {
            tmpTotalSum += tmpAdderGroup.getPressureZAdder().getSum();
        }
        return tmpTotalSum;
    }
    
    /**
     * Return accumulated (total) sum of all adders for G value of PNHLN 
     * integration
     * NOTE: NO checks are performed.
     * 
     * @return Return accumulated (total) sum of all adders for G value of PNHLN 
     * integration
     */
    public double getAccumulatedPnhln_G_AddersSum() {
        double tmpTotalSum = 0.0;
        for (AdderGroup tmpAdderGroup : this.adderGroups) {
            tmpTotalSum += tmpAdderGroup.getPnhln_G_Adder().getSum();
        }
        return tmpTotalSum;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get/set)">
    /**
     * Executor service
     * 
     * @return Executor service
     */
    public ExecutorService getExecutorService() {
        return this.executorService;
    }

    /**
     * Executor service
     * NOTE: NO checks are performed.
     * 
     * @param anExecutorService Executer service
     */
    public void setExecutorService(ExecutorService anExecutorService) {
        this.executorService =  anExecutorService;
    }
    
    /**
     * Adder groups
     * 
     * @return Adder groups
     */
    public AdderGroup[] getAdderGroups() {
        return this.adderGroups;
    }
    
    /**
     * Adder groups
     * NOTE: NO checks are performed.
     * 
     * @param anAdderGroups Adder groups
     */
    public void setAdderGroups(AdderGroup[] anAdderGroups) {
        this.adderGroups = anAdderGroups;
    }
    // </editor-fold>
    
}
