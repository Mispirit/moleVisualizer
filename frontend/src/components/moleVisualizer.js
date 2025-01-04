import React, { useState, useRef } from "react";
import axios from "axios";
import * as $3Dmol from "3dmol"; // Import 3Dmol.js for 3D visualization
import Head from "./header";
import "../App.css";
import Footer from "./footer";
import { IconButton, CircularProgress } from "@mui/material";
import { Search, Clear } from "@mui/icons-material";

function MoleVisualizer() {
  const [smiles, setSmiles] = useState("");
  const [properties, setProperties] = useState(null);
  const [error, setError] = useState("");
  const [mol2D, setMol2D] = useState("");
  const [compoundData, setCompoundData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [isSearching, setIsSearching] = useState(false); // Toggle between Search and Clear button

  const mol3DRef = useRef(null); // Reference to the 3Dmol viewer container
  const [viewer, setViewer] = useState(null); // To store the viewer instance

  // Fetch compound name
  const fetchCompoundName = async () => {
    setLoading(true); // Set loading to true when fetching data
    try {
      const response = await fetch(
        `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${smiles}/property/IUPACName/JSON`
      );
      if (!response.ok) {
        throw new Error("Failed to fetch compound data.");
      }
      const data = await response.json();
      const name = data?.PropertyTable?.Properties?.[0]?.IUPACName;
      if (!name) {
        throw new Error("Invalid compound data structure.");
      }
      setCompoundData({ name });
    } catch (error) {
      console.error("Error fetching compound name:", error);
      alert("Failed to fetch compound name!");
    } finally {
      setLoading(false); // Set loading to false once data is fetched or fails
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    
    // Prevent submitting if input is empty
    if (!smiles) return;

    setIsSearching(true); // Change button to "Clear" after searching
    try {
      const response = await axios.post("http://localhost:5001/api/calculate", {
        smiles,
      });
      setProperties(response.data);
      setError("");

      // Fetch compound name
      await fetchCompoundName();

      // Fetch 2D structure image
      const image2D = await axios.post(
        "http://localhost:5001/api/render2d",
        { smiles },
        { responseType: "blob" }
      );
      setMol2D(URL.createObjectURL(image2D.data));

      // Render 3D structure
      render3D(smiles);
    } catch (err) {
      if (err.response) {
        setError("Error calculating properties: " + err.response.data.error);
      } else if (err.request) {
        setError(
          "No response from server. Please check if the backend is running."
        );
      } else {
        setError("An unexpected error occurred: " + err.message);
      }
      setProperties(null);
    }
  };

  // Render 3D molecular structure using 3Dmol.js
  const render3D = (smiles) => {
    if (mol3DRef.current) {
      const viewerInstance = $3Dmol.createViewer(mol3DRef.current, {
        backgroundColor: "white",
      });

      // Store the viewer instance
      setViewer(viewerInstance);

      axios
        .get(
          `https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(
            smiles
          )}/sdf`
        )
        .then((res) => {
          if (!res.data) {
            throw new Error("Empty SDF data received.");
          }
          viewerInstance.addModel(res.data, "sdf");
          viewerInstance.setStyle({}, { stick: { radius: 0.2 }, sphere: { scale: 0.3 } });
          viewerInstance.zoomTo();
          viewerInstance.render();
        })
        .catch((err) => {
          console.error("Error fetching 3D structure:", err);
          alert("Failed to render 3D structure.");
        });
    }
  };

  // Clear function to reset the state
  const clear = () => {
    setSmiles(""); // Clear the input
    setProperties(null); // Clear properties
    setMol2D(""); // Clear 2D structure
    setCompoundData(null); // Clear compound name
    setError(""); // Clear error message
    setIsSearching(false); // Reset to show the Search button
    if (viewer) {
      viewer.clear(); // Clear the 3Dmol viewer
    }
  };

  return (
    <div className="App">
      <Head />
      <form onSubmit={handleSubmit} className="max-w-[1200px] extrahead">
        <div className="max-w-[90%] flex justify-between lg:max-w-[800px] extrahead p-[15px] bg-white rounded 3xl shadow-lg">
          <input
            type="text"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            placeholder="Enter SMILE STRING e.g., C1CCCCC1"
            className="w-[90%] text-black focus:outline-none focus:ring-0"
            disabled={isSearching} // Disable input during search
          />
          <IconButton
            onClick={isSearching ? clear : handleSubmit}
            disabled={isSearching && !smiles} // Disable button if searching is not in progress
            style={{
              backgroundColor: "#007bff", 
              color: "white", 
              marginLeft: "10px",
            }}
          >
            {isSearching ? <Clear /> : <Search />}
          </IconButton>
        </div>
      </form>

      <br />

      {error && <p style={{ color: "red" }}>{error}</p>}

      {properties && (
        <div className="bg-white">
          <div className="bg-white block lg:flex lg:justify-between max-w-[1200px] p-[20px] extrahead">
            <div className="bg-white block lg:flex lg:justify-between w-[100%] ">
              <div className="w-[90%] lg:w-[46%] shadow-lg extrahead lg:pl-[20px]">
                <div className="w-[100%] overflow-hidden ">
                  <p className="text-3xl font-bold">3D Structure:</p>
                  <div
                    ref={mol3DRef}
                    className="relative top-[0] w-full min-w-[300px] min-h-[300px] lg:min-w-[400px] lg:min-h-[400px] overflow-hidden z-index-[-1000]"
                  ></div>
                </div>
                <div className="w-[100%]">
                  <p className="text-3xl font-bold">2D Structure:</p>
                  {mol2D && (
                    <img
                      src={mol2D}
                      alt="2D molecular structure"
                      className="extrahead"
                    />
                  )}
                </div>
              </div>

              <div className="w-[90%] lg:w-[49%] shadow-lg extrahead p-[20px] bg-white">
                {compoundData && compoundData.name && (
                  <h1 className="text-3xl capitalize">{compoundData.name}</h1>
                )}
                <br />
                <p className="text-2xl font-bold pb-[10px]">Properties:</p>
                <ul className="listyle">
                  <li>
                    <strong>Molecular Weight:</strong>{" "}
                    {properties.molecular_weight}
                  </li>
                  <li>
                    <strong>LogP:</strong> {properties.logP}
                  </li>
                  <li>
                    <strong>TPSA:</strong> {properties.TPSA}
                  </li>
                  <li>
                    <strong>Hydrogen Bond Acceptors:</strong>{" "}
                    {properties.num_h_acceptors}
                  </li>
                  <li>
                    <strong>Hydrogen Bond Donors:</strong>{" "}
                    {properties.num_h_donors}
                  </li>
                  <li>
                    <strong>Rotatable Bonds:</strong>{" "}
                    {properties.num_rotatable_bonds}
                  </li>
                  <li>
                    <strong>Molecular Formula:</strong> {properties.formula}
                  </li>
                </ul>
              </div>
            </div>
          </div>
        </div>
      )}

      <br />
      {loading && <CircularProgress />}
      <div className="fixed bottom-[0] w-[100%]"><Footer /></div>
    </div>
  );
}

export default MoleVisualizer;
