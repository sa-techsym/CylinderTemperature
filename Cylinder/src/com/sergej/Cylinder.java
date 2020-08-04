package com.sergej;

public class Cylinder {
	private double _inner_radius = .05, _outer_radius = .06, _length = 1.0, _heat_conduct = 70.0, _temp_conduct = 2.3e-5;
		
	public float inner_radius() {
		return (float) _inner_radius;
		}
		
	public float outer_radius() {
		return (float) _outer_radius;
		}
		
	public float length() {
		return (float) _length;
		}
		
	public float heat_conduct() {
		return (float) _heat_conduct;
		}
			
	public float temp_conduct() {
		return (float) _temp_conduct;
		}
	/*_heat_transfer_flow, _heat_transfer_envr*/;
	}
