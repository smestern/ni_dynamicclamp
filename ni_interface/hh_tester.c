// A hodgkin huxley neuron model tester
// used for testing the neuron interface by providing a HH model that spikes independently

// adapted from the brilliant work of swharden: https://github.com/swharden/HHSharp

#include <math.h>
#include <stdio.h>

// Function pointer type for updating time constants
typedef void (*UpdateTimeConstantsFunc)(void* gate, double Vm);

// Base voltage gate structure
typedef struct {
    double alpha;  // rate of opening (units = fraction / millisecond)
    double beta;   // rate of closing (units = fraction / millisecond)
    double activation;  // activation fraction (range 0-1)
    UpdateTimeConstantsFunc updateTimeConstants;
} VoltageGate;

// Helper functions for voltage gates
double get_activation_change_per_ms(VoltageGate* gate) {
    return gate->alpha * (1 - gate->activation) - gate->beta * gate->activation;
}

void set_infinite_state(VoltageGate* gate) {
    gate->activation = gate->alpha / (gate->alpha + gate->beta);
}

void step_forward(VoltageGate* gate, double stepSizeMs) {
    gate->activation += get_activation_change_per_ms(gate) * stepSizeMs;
}

// Specific gate update functions
void vgsc_activation_update(void* gate, double Vm) {
    VoltageGate* g = (VoltageGate*)gate;
    g->alpha = 0.1 * ((25 - Vm) / (exp((25 - Vm) / 10) - 1));
    g->beta = 4 * exp(-Vm / 18);
}

void vgsc_inactivation_update(void* gate, double Vm) {
    VoltageGate* g = (VoltageGate*)gate;
    g->alpha = 0.07 * exp(-Vm / 20);
    g->beta = 1 / (exp((30 - Vm) / 10) + 1);
}

void vgkc_activation_update(void* gate, double Vm) {
    VoltageGate* g = (VoltageGate*)gate;
    g->alpha = 0.01 * ((10 - Vm) / (exp((10 - Vm) / 10) - 1));
    g->beta = 0.125 * exp(-Vm / 80);
}

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

// Initialize HH Model
void hh_model_init(HHModel* model, double initialVoltageOffset) {
    model->ENa = 115;
    model->EK = -12;
    model->EKleak = 10.6;
    model->gNa = 120;
    model->gK = 36;
    model->gKleak = 0.3;
    model->Cm = 1;
    model->Vm = initialVoltageOffset;
    
    // Setup gates with their respective update functions
    model->m.updateTimeConstants = vgsc_activation_update;
    model->h.updateTimeConstants = vgsc_inactivation_update;
    model->n.updateTimeConstants = vgkc_activation_update;
    
    // Initialize gates
    model->m.updateTimeConstants(&model->m, model->Vm);
    model->h.updateTimeConstants(&model->h, model->Vm);
    model->n.updateTimeConstants(&model->n, model->Vm);
    
    set_infinite_state(&model->m);
    set_infinite_state(&model->h);
    set_infinite_state(&model->n);
}

// Update all gate time constants
void update_all_gate_time_constants(HHModel* model) {
    model->m.updateTimeConstants(&model->m, model->Vm);
    model->h.updateTimeConstants(&model->h, model->Vm);
    model->n.updateTimeConstants(&model->n, model->Vm);
}

// Update currents and voltage
void update_currents_and_voltage(HHModel* model, double stimulusCurrent, double deltaTms) {
    model->INa = pow(model->m.activation, 3) * model->gNa * model->h.activation * (model->Vm - model->ENa);
    model->IK = pow(model->n.activation, 4) * model->gK * (model->Vm - model->EK);
    model->IKleak = model->gKleak * (model->Vm - model->EKleak);
    model->Isum = stimulusCurrent - model->INa - model->IK - model->IKleak;
    model->Vm += deltaTms * model->Isum / model->Cm;
}

// Step the model forward in time
void hh_model_step(HHModel* model, double stimulusCurrent, double stepSizeMs) {
    update_all_gate_time_constants(model);
    update_currents_and_voltage(model, stimulusCurrent, stepSizeMs);
    step_forward(&model->m, stepSizeMs);
    step_forward(&model->h, stepSizeMs);
    step_forward(&model->n, stepSizeMs);
}

