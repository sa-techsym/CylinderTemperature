package com.sergej;

import java.util.ArrayList;
import java.util.function.Function;
import cern.jet.math.Bessel;

public class CylinderTemperature {
	//private double _inner_radius, _outer_radius, _length, _heat_conduct, _temp_conduct; /*_heat_transfer_flow, _heat_transfer_envr*/;
	Cylinder _cylinder;
	
	CylinderTemperature(Cylinder source) {
		_cylinder = source;
		}
	
	ArrayList <Double> _bessel_characteristic, _trigonometric_characteristic;
	
	public CylinderTemperature calculate(double flow_heat_transfer, double envr_heat_transfer) {
		_flow_heat_transfer = flow_heat_transfer;
		_envr_heat_transfer = envr_heat_transfer;
		
		_bessel_characteristic 
				= solve(x ->
						  (x * _cylinder.heat_conduct() * Bessel.j1(x * _cylinder.outer_radius()) - _envr_heat_transfer 
								  * Bessel.j0(x * _cylinder.outer_radius())
								  ) 
						* (x * _cylinder.heat_conduct() * Bessel.y1(x * _cylinder.inner_radius()) + _flow_heat_transfer 
								  * Bessel.y0(x * _cylinder.inner_radius())
								  )	- 
						  (x * _cylinder.heat_conduct() * Bessel.y1(x * _cylinder.outer_radius()) - _envr_heat_transfer 
								  * Bessel.y0(x * _cylinder.outer_radius())
								  ) 
						* (x * _cylinder.heat_conduct() * Bessel.j1(x * _cylinder.inner_radius()) + _flow_heat_transfer 
								  * Bessel.j0(x * _cylinder.inner_radius())
								  )
						);
		for (double item : _bessel_characteristic)
			System.out.println(item);
		
		_trigonometric_characteristic 
				= solve(x -> {
						double phase = Math.atan(x * _cylinder.heat_conduct() / _envr_heat_transfer);
						return  x * _cylinder.heat_conduct() * Math.cos(x * _cylinder.length() + phase) + _envr_heat_transfer 
								* Math.sin(x * _cylinder.length() + phase); 
								}
						); 
		for (double item : _trigonometric_characteristic)
			System.out.println(item);
	
		return this;
		}
	
	private double value(double t, double r, double x) {
		double result = 20.0; 
		System.out.print("->");
		for (double bessel_item : _bessel_characteristic)
			for (double trigonometric_item : _trigonometric_characteristic)
				result += v(t, bessel_item, trigonometric_item) * p(r, bessel_item) * s(x, trigonometric_item)
						/ (r * n(trigonometric_item) * z(bessel_item));
			return result;
		}
	
	private double _flow_heat_transfer, _envr_heat_transfer;
	
	private double n(double /*trigonometric*/ characteristic) {
		double phase = Math.atan(characteristic * _cylinder.heat_conduct() / _envr_heat_transfer);
		return 
				.5 * (_cylinder.length() - (Math.sin(characteristic * _cylinder.length() + phase) * Math.cos(characteristic 
						* _cylinder.length() + phase) - Math.sin(phase) * Math.cos(phase)) / characteristic);
		}
	
	private double z(double /*bessel*/ characteristic) {
		Function <Double, Double> value 
				= x -> 
						(Bessel.j0(x) + a(characteristic) * Bessel.y0(x)) * (Bessel.j0(x) + a(characteristic) 
								* Bessel.y0(x))
					  + (Bessel.j1(x) + a(characteristic) * Bessel.y1(x)) * (Bessel.j1(x) + a(characteristic) 
					  			* Bessel.y1(x)); 
		return 
				.5 * _cylinder.outer_radius() * _cylinder.outer_radius() * value.apply(characteristic * _cylinder.outer_radius()) 
			  - .5 * _cylinder.inner_radius() * _cylinder.inner_radius() * value.apply(characteristic * _cylinder.inner_radius());
		}

	private double a(double /*bessel*/ characteristic) {
		double freq = characteristic * _cylinder.inner_radius();
		return 
				-1.0 * (characteristic * _cylinder.heat_conduct() * Bessel.j1(freq) + _flow_heat_transfer * Bessel.j0(freq)) 
				     / (characteristic * _cylinder.heat_conduct() * Bessel.y1(freq) + _flow_heat_transfer * Bessel.y0(freq)); 
		}
	
	private double p(double r, double /*bessel*/ characteristic) {
		return r * (Bessel.j0(characteristic * r) + a(characteristic) * Bessel.y0(characteristic * r));
		}
	
	private double s(double x, double /*trigonometric*/ characteristic) {
		double phase = Math.atan(characteristic * _cylinder.heat_conduct() / _envr_heat_transfer);
		return Math.sin(characteristic * x + phase);
		}
	
	private double v(double t, double bessel_characteristic, double trigonometric_characteristic) {
		/* s(x) function integral from 0 to full temperature profile length */
		double 
			phase = Math.atan(trigonometric_characteristic * _cylinder.heat_conduct() / _envr_heat_transfer),
			trigonometric_integral = (Math.cos(phase) - Math.cos(trigonometric_characteristic * _cylinder.length() + phase)) 
					/ trigonometric_characteristic;
		/*System.out.print(trigonometric_integral);
		 *System.out.print(" ");*/
		
		/* r * p(r) function integral from inner radius to outer radius 
		 * bfjy(r0*hp[i],ji0,ji1,yi0,yi1); bfjy(r1*hp[i],j00,j01,y00,y01); return((r1*j01-r0*ji1+a[i]*(r1*y01-r0*yi1))/hp[i]);*/
		double bessel_integral = 
			    (_cylinder.outer_radius() * Bessel.j1(_cylinder.outer_radius() * bessel_characteristic) - _cylinder.inner_radius() 
						* Bessel.j1(_cylinder.inner_radius() * bessel_characteristic) 
			  + a(bessel_characteristic) 
			  * (_cylinder.outer_radius() * Bessel.y1(_cylinder.outer_radius() * bessel_characteristic) - _cylinder.inner_radius() 
						* Bessel.y1(_cylinder.inner_radius() * bessel_characteristic))
			    		) / bessel_characteristic;
		/*System.out.print(bessel_integral);
		 *System.out.print(" ");*/
		
		double 	sum = bessel_characteristic * bessel_characteristic + trigonometric_characteristic * trigonometric_characteristic;
	 	/* v[j][i]=fta(tay)*is(xl)*(t0-tc0)*ic()+((tcv-tc0)*al1*fr(r0)+(tcf-tc0)*al2*fr(r1))*is(xl)*(1.0-fta(tay))/(la*h[j][i]); */
		
		Function <Double, Double> u = x -> { if (x > 60 ) return 0.0; else return Math.exp(-1 * x); };
				
	 	double debug_a = u.apply(_cylinder.temp_conduct() * sum * t) * trigonometric_integral * (/*t0*/20 - /*tc0*/20) * bessel_integral 
	 			+ (
	 					(/*tcv*/50 - /*tc0*/20) * _flow_heat_transfer * p(_cylinder.inner_radius(), bessel_characteristic) 
	 				  + (/*tcf*/20 -/*tc0*/ 20) * _envr_heat_transfer * p(_cylinder.outer_radius(), bessel_characteristic)
	 				  ) 
	 			* trigonometric_integral * (1.0 - u.apply(_cylinder.temp_conduct() * sum * t)) / (_cylinder.heat_conduct() * sum);
	 	/*System.out.print(debug_a);
		 *System.out.println(" ");*/
		
	 	return debug_a;
		}

	private ArrayList <Double> solve(Function <Double, Double> equation) {
		ArrayList <Double> characteristic = new ArrayList <Double> ();
		
		double current = 0.000001f;
		
		for (int i = 0; i < 10; i++) {
			double a = equation.apply(current);
			while (Math.signum(equation.apply(current)) * Math.signum(equation.apply(current += 0.05f)) == 1.0f)
				;
			double left = current - 0.05f, right = current;
		
			while (right - left > 0.01) {
				double temp = (right + left) / 2;
				if (equation.apply(right) * equation.apply(temp) < 0)
					left = temp;
				else
					right = temp;
				}
			characteristic.add((right + left) / 2);
			}
		
		return characteristic;
		}
	
	public static void main(String[] args) {
		CylinderTemperature cylinder = new CylinderTemperature(new Cylinder()).calculate(2000,  20);
		for (int i = 0; i < 20000; i += 10)
			System.out.println(cylinder.value(i, .06, 0.9));
		}
	}


