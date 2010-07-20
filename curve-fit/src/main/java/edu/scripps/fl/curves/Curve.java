/*
 * Copyright 2010 The Scripps Research Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package edu.scripps.fl.curves;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang.builder.ToStringBuilder;

/**
 * 
 * @author Mark Southern (southern at scripps dot edu)
 * 
 */
public class Curve {
	private List<Double>	concentrations	= new ArrayList();
	private Double			curveClass;
	private String			curveDescription;
	private Double			EC50;
	private Double			hillSlope;
	private Double			IC50;
	private Double			logEC50;
	private List<Boolean>	mask			= new ArrayList();
	private Boolean			masked;
	private Double			maxResponse;
	private Double			pHill;
	private Double			R2;
	private Double			rank;
	private Double			responseMax;
	private Double			responseMin;
	private Double			responseRange;
	private List<Double>	responses		= new ArrayList();
	private String			signalDirection;
	private Double			SYX;
	private Double			YInflection;
	private Double			YZero;

	public void add(Double response, Double concentration) {
		responses.add(response);
		concentrations.add(concentration);
	}

	public List<Double> getConcentrations() {
		return concentrations;
	}

	public Double getCurveClass() {
		return curveClass;
	}

	public String getCurveDescription() {
		return curveDescription;
	}

	public Double getEC50() {
		return EC50;
	}

	public Double getHillSlope() {
		return hillSlope;
	}

	public Double getIC50() {
		return IC50;
	}

	public Double getLogEC50() {
		return logEC50;
	}

	public List<Boolean> getMask() {
		return mask;
	}

	public Boolean getMasked() {
		return masked;
	}

	public Double getMaxResponse() {
		return maxResponse;
	}

	public Double getPHill() {
		return pHill;
	}

	public Double getR2() {
		return R2;
	}

	public Double getRank() {
		return rank;
	}

	public Double getResponseMax() {
		return responseMax;
	}

	public Double getResponseMin() {
		return responseMin;
	}

	public Double getResponseRange() {
		return responseRange;
	}

	public List<Double> getResponses() {
		return responses;
	}

	public String getSignalDirection() {
		return signalDirection;
	}

	public Double getSYX() {
		return SYX;
	}

	public Double getYInflection() {
		return YInflection;
	}

	public Double getYZero() {
		return YZero;
	}

	public Boolean isMasked() {
		return masked;
	}

	public void setConcentrations(List<Double> concentrations) {
		this.concentrations = concentrations;
	}

	public void setCurveClass(Double curveClass) {
		this.curveClass = curveClass;
	}

	public void setCurveDescription(String curveDescription) {
		this.curveDescription = curveDescription;
	}

	public void setEC50(Double ec50) {
		EC50 = ec50;
	}

	public void setHillSlope(Double hillSlope) {
		this.hillSlope = hillSlope;
	}

	public void setIC50(Double ic50) {
		IC50 = ic50;
	}

	public void setLogEC50(Double logEC50) {
		this.logEC50 = logEC50;
	}

	public void setMask(List<Boolean> mask) {
		this.mask = mask;
	}

	public void setMasked(Boolean masked) {
		this.masked = masked;
	}

	public void setMaxResponse(Double maxResponse) {
		this.maxResponse = maxResponse;
	}

	public void setPHill(Double hill) {
		pHill = hill;
	}

	public void setR2(Double r2) {
		R2 = r2;
	}

	public void setRank(Double rank) {
		this.rank = rank;
	}

	public void setResponseMax(Double responseMax) {
		this.responseMax = responseMax;
	}

	public void setResponseMin(Double responseMin) {
		this.responseMin = responseMin;
	}

	public void setResponseRange(Double responseRange) {
		this.responseRange = responseRange;
	}

	public void setResponses(List<Double> responses) {
		this.responses = responses;
	}

	public void setSignalDirection(String signalDirection) {
		this.signalDirection = signalDirection;
	}

	public void setSYX(Double syx) {
		SYX = syx;
	}

	public void setYInflection(Double inflection) {
		YInflection = inflection;
	}

	public void setYZero(Double zero) {
		YZero = zero;
	}

	public String toString() {
		return ToStringBuilder.reflectionToString(this);
	}
}