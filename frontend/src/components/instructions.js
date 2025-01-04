import React, { useState } from "react";
import Mollogo from "../images/mollogo.png";
import "../styles/extrahead.css";
import { Link } from "react-router-dom";
import Footer from "./footer";
import "../App.css"
function Instructions() {
  return (
    <div className=" ">
              <div className="bg-black w-[100%]">

      <div className="flex justify-between p-[20px] max-w-[1000px] extrahead bg-black">
        {" "}
        <img src={Mollogo} className="w-[180px] " />
        <Link to="/">
          {" "}
          <p className="text-accent">Visualizer</p>
        </Link>
      </div></div>
      <br />
      <div className="instructions-container w-[90%] max-w-[1000px] extrahead ">
        <h2 className="text-2xl font-bold text-accent-hover h2font">
          Instructions for Using the Molecule Visualizer
        </h2>
        <p className="instructionwlc">
          Welcome to my Molecule Visualizer Project! This tool allows you to
          explore the 2D and 3D structures of molecules by entering their SMILES
          (Simplified Molecular Input Line Entry System) string. Follow these
          simple steps to get started:
        </p>

        <br/>

        <ol className="listed">
          <li>
            <strong>Enter the SMILES String:</strong>
            <ul>
              <li>
                A SMILES string is a text representation of a chemical
                structure. Example: <code>C1CCCCC1</code> for Cyclohexane.
              </li>
              <li>Type or paste a valid SMILES string into the input field.</li>
            </ul>
          </li>

          <li>
            <strong>Click the 'Search' Button:</strong>
            <ul>
              <li>
                After entering a SMILES string, click the{" "}
                <strong>Search</strong> Icon button to load the molecular
                properties, 2D image, and 3D visualization.
              </li>
              <li>The app will fetch relevant information and display it.</li>
            </ul>
          </li>

          <li>
            <strong>View the Molecular Structure:</strong>
            <ul>
              <li>
                <strong>2D Structure:</strong> The 2D structure of the molecule
                will be displayed below the 3D view.
              </li>
              <li>
                <strong>3D Structure:</strong> A 3D visualization of the
                molecule will appear below the input field.
              </li>
            </ul>
          </li>

          <li>
            <strong>View the Properties of the Molecule:</strong>
            <ul>
              <li>
                Molecular properties such as Molecular Weight, LogP, TPSA,
                Hydrogen Bond Acceptors/Donors, and more will be listed below
                the 3D structure.
              </li>
              <li>The IUPAC name of the compound will also be displayed.</li>
            </ul>
          </li>

          <li>
            <strong>Clear the Inputs:</strong>
            <ul>
              <li>
                If you want to enter a new SMILES string or reset the
                application, click the <strong>Clear (X)</strong> button. This
                will reset the input field, remove the displayed properties, and
                clear both the 2D and 3D structures.
              </li>
            </ul>
          </li>

          <li>
            <strong>Error Handling:</strong>
            <ul>
              <li>
                If an error occurs, an error message will be displayed. Make
                sure the SMILES string is valid and try again.
              </li>
            </ul>
          </li>
        </ol>

        <hr />

        <h2 className="text-2xl font-bold text-accent-hover h2font">
        Extra Tips:</h2>
        <ul className="listed">
          <li>
            <strong>Valid SMILES Strings:</strong> Ensure the SMILES string you
            input is correctly formatted. If you're unsure about the validity,
            use an online SMILES generator or checker.
          </li>
          <li>
            <strong>Zoom Model:</strong> You can interact with the 3D structure
            using your mouse. Zoom in or out by scrolling.
          </li>
          <li>
            <strong>Slow Response:</strong> In case the application is slow to
            respond, it might be due to high traffic or large molecule size.
            Please be patient while it processes the data.
          </li>
        </ul>
      </div>


      <br/>
      <br/>
      <br/>
      <div className="fixed bottom-[0] w-[100%]">
        <Footer />
      </div>
    </div>
  );
}

export default Instructions;
