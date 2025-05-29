# Motor Design Tool (MDT)

A comprehensive software tool for designing and simulating electric motors based on end-use case requirements, specifically tailored for electric vehicle applications.

## Overview

The Motor Design Tool enables users to design and simulate electric motors by starting with vehicle parameters and performance specifications. The tool derives electric motor specifications from vehicle requirements and provides both analytical solutions and Finite Element Method (FEMM) analysis for iterative design refinement.

## Architecture

### Frontend
- **Technology**: Flutter Web Stack
- **Hosting**: AWS S3 bucket with CloudFront distribution
- **Purpose**: Web-based user interface for motor design and simulation

### API Server
- **File**: `server_multi_db.py`
- **Function**: REST API endpoint connecting frontend and backend
- **Features**: 
  - Multiprocessing for background simulation tasks
  - Result polling capability
  - MongoDB database integration

### Backend

#### 1. Vehicle Dynamics Module
Derives motor requirements from vehicle specifications.

**Files:**
- `vehicle_dynamics.py` - Contains VehicleDynamics class
- `vehicle_dynamics_methods.py` - Supporting functions

**Outputs:**
- Continuous mode and plot
- Peak mode and plot  
- Voltage specifications

#### 2. Motor Wiz Module
Calculates motor geometry, winding patterns, parameters, cost, and operating points.

**Components:**

**Materials Subsystem:**
- `magnet.ini` - Magnet material parameters
- `wire.ini` - Wire specifications
- `steel.ini` - Steel material properties
- `htc.ini` - Thermal material properties
- `load_materials.py` - Material loading functions
- `thermal.py` - Thermal data management

**Topology Subsystem:**
- `inner_rotor.py` - Inner rotor calculations for SPMSM and IPMSM motors
- `motorwiz.py` - Main MotorWiz class and orchestration

**Outputs:**
- Motor DXF files
- Winding DXF files
- Internal calculated variables
- Performance specifications
- Material cost and weight analysis

#### 3. Test Design Module
Converts motor geometry into DXF files for FEMM analysis.

**Core Files:**
- `dxf.py` - Motor cross-section and winding diagram generation
- `run_femm.py` - FEMM simulation integration

**Topology Support:**
- `IPMSM_dxf.py` - Interior Permanent Magnet Synchronous Motor drawing
- `SPMSM_dxf.py` - Surface Permanent Magnet Synchronous Motor drawing
- `IPMSM_femm.py` - IPMSM simulation
- `SPMSM_femm.py` - SPMSM simulation

**FEMM Analysis Outputs:**
- Air gap flux plots and values (Bg_avg)
- Magnetic flux density plots
- Inductance parameters (Ld, Lq)
- Magnetic flux linkage (psi_m)

## Supported Motor Types

- **IPMSM** (Interior Permanent Magnet Synchronous Motor)
- **SPMSM** (Surface Permanent Magnet Synchronous Motor)
- **Radial topology** configuration

## Workflow

1. **Input vehicle parameters** and performance specifications
2. **Generate motor requirements** using Vehicle Dynamics module
3. **Design motor geometry** with Motor Wiz calculations
4. **Create DXF files** for motor and winding patterns
5. **Run FEMM simulations** for performance validation
6. **Iterate design** based on simulation results

## Dependencies

- Python backend with multiprocessing support
- FEMM (Finite Element Method Magnetics)
- MongoDB database
- Flutter Web framework
- AWS services (S3, CloudFront)

## Getting Started

The tool is designed for iterative use, allowing engineers to refine motor designs through multiple simulation cycles until achieving optimal performance for their specific electric vehicle application.
