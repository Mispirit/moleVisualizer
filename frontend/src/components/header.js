?import React, { useState } from "react";
import axios from "axios";
import * as $3Dmol from "3dmol"; // Import 3Dmol.js for 3D visualization
import Mollogo from "../images/mollogo.png";
import "../styles/extrahead.css";
import { Link } from "react-router-dom";

function head() {
  return (
    <div className="bg-black  p-[20px] h-[300px] ">
      <div className="flex justify-between max-w-[1000px] extrahead">
        {" "}
        <img src={Mollogo} className="w-[120px] h-[20px] lg:w-[180px] lg:h-[30px]" />
       <Link to="/Instructions"> <p className="text-accent">Instructions</p></Link>
      </div>

      <br/>
      <br/>
      <br/>

      {/* About Web*/}
<h1 className="text-white text-center text-3xl lg:text-4xl font-bold leading-[2.5rem]"> 2D/3D Chemical Structure <br/> Visualizer and Properties Tool</h1>
   <p> Please hold on for 30 Secs</p> </div>
  );
}

export default head;
