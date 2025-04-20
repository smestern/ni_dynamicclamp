#ifndef HH_TESTER_H
#define HH_TESTER_H

#ifdef __cplusplus
extern "C" {
#endif

// Function pointer type for updating time constants
typedef void (*UpdateTimeConstantsFunc)(void* gate, double Vm);

// Base voltage gate structure
typedef struct {
    double alpha;  // rate of opening (units = fraction / millisecond)
    double beta;   // rate of closing (units = fraction / millisecond)
    double activation;  // activation fraction (range 0-1)
    UpdateTimeConstantsFunc updateTimeConstants;
} VoltageGate;

// HH Model structure
typedef struct {
    double ENa;    // sodium reversal potential (mV)
    double EK;     // potassium reversal potential (mV)
    double EKleak; // leak reversal potential (mV)
    double gNa;    // sodium conductance (mS/cm²)
    double gK;     // potassium conductance (mS/cm²)
    double gKleak; // leak conductance (mS/cm²)
    double Cm;     // membrane capacitance (µF/cm²)
    
    VoltageGate m; // sodium activation gate
    VoltageGate h; // sodium inactivation gate
    VoltageGate n; // potassium activation gate
    
    double INa;    // sodium current (µA/cm²)
    double IK;     // potassium current (µA/cm²)
    double IKleak; // leak current (µA/cm²)
    double Isum;   // total current (µA/cm²)
    double Vm;     // membrane voltage (mV)
} HHModel;

// Helper functions for voltage gates
double get_activation_change_per_ms(VoltageGate* gate);
void set_infinite_state(VoltageGate* gate);
void step_forward(VoltageGate* gate, double stepSizeMs);

// Gate update functions
void vgsc_activation_update(void* gate, double Vm);
void vgsc_inactivation_update(void* gate, double Vm);
void vgkc_activation_update(void* gate, double Vm);

// HH Model functions
void hh_model_init(HHModel* model, double initialVoltageOffset);
void update_all_gate_time_constants(HHModel* model);
void update_currents_and_voltage(HHModel* model, double stimulusCurrent, double deltaTms);
void hh_model_step(HHModel* model, double stimulusCurrent, double stepSizeMs);

#ifdef __cplusplus
}
#endif

#endif // HH_TESTER_H